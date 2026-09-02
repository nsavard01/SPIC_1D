#include <vector>
#include "particles/charged_particle.hpp"
#include "globals/mpi_vars.hpp"
#include "globals/constants.hpp"
#include <omp.h>
#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------------------------
// Implicit CIC (cloud-in-cell) mover.
//
// The cells a particle travels through are bounded by the field nodes rather than by the
// potential nodes, which keeps the cell index a continuous function of position.  In
// logical coordinates the field nodes sit on the integers 0..number_cells and the
// potential/charge nodes sit half way between them, at cell centre c + 1/2.  A particle
// therefore sees a field that varies linearly across its cell,
//
//     a(s) = [ dphi_left * (1 - s) + dphi_right * s ] * (q/m) / dx_cell ,   s = xi - cell,
//
// working from the potential difference across each cell edge rather than the field itself,
// so that the cell width used to reach logical coordinates cancels against the length each
// edge value spans.  On a uniform grid the two are the same; on a stretched grid only this
// form keeps the work done on the particles equal to the change in field energy.
//
// Charge is deposited with the quadratic B-spline that is the integral of that linear
// weighting, which is what makes the scheme energy conserving.
// ---------------------------------------------------------------------------------------

namespace {

// convergence criterion of the Picard iteration, in logical (cell) units
constexpr double picard_tolerance = 1e-10;
// safety valve, the iteration converges in a handful of passes in normal operation
constexpr int max_picard_iterations = 50;
// safety valve on the number of cell crossings within one time step
constexpr int max_sub_steps = 1000;

// Deposit one particle onto the cell-centred charge nodes with the quadratic B-spline.
// "cell" is the cell the particle is tracked in, which for a particle sitting exactly on
// the periodic seam is not int(xi_p), so it is passed in explicitly.
inline void deposit_quadratic(std::vector<double>& work_space, const double xi_p, const int cell,
    const int left_boundary, const int right_boundary, const int number_cells) {
    const double d = xi_p - double(cell);
    const double weight_left = 0.5 * (1.0 - d) * (1.0 - d);
    const double weight_right = 0.5 * d * d;
    work_space[cell] += 0.5 + d - d * d;
    if (cell > 0) {
        work_space[cell-1] += weight_left;
    } else {
        switch (left_boundary) {
            case 1:
            case 4:
                // Dirichlet, the part of the cloud past the wall is subtracted so the deposit
                // stays consistent with the half-cell gradient the field solver uses there
                work_space[0] -= weight_left;
                break;
            case 2:
                // reflective/symmetric, the part past the wall folds back in
                work_space[0] += weight_left;
                break;
            case 3:
                work_space[number_cells-1] += weight_left;
                break;
        }
    }
    if (cell < number_cells-1) {
        work_space[cell+1] += weight_right;
    } else {
        switch (right_boundary) {
            case 1:
            case 4:
                work_space[number_cells-1] -= weight_right;
                break;
            case 2:
                work_space[number_cells-1] += weight_right;
                break;
            case 3:
                work_space[0] += weight_right;
                break;
        }
    }
}

// One sub-step of a particle inside a single cell, in the local coordinate s in [0,1].
// Returns true when the particle reaches a cell edge within del_tau, in which case del_tau is
// replaced by the time actually taken, s_f is set exactly onto the edge, and local_boundary
// reports which edge was reached (0 = left, 1 = right).
//
// Whether an edge is reached is decided from the energy available at the start of the
// sub-step, not from where the trajectory happens to end.  An endpoint test cannot see a
// particle that leaves the cell and returns within the same sub-step, and would then let the
// iteration below extrapolate this cell's field beyond its own boundary.
inline bool CIC_sub_step(const double s_i, const double v_i, double& del_tau, double& s_f, double& v_f,
    const CIC_cell_coefficients& cell, int& local_boundary) {

    const double inv_dx = cell.inv_dx;
    const double b_left = cell.q_dphi_left;
    const double b_right = cell.q_dphi_right;
    const double gradient = b_right - b_left;
    const double v_i_sqr = v_i * v_i;

    // Work done reaching s, in the units where v^2(s) = v_i^2 + work_to(s).  The cell width
    // cancels because q_dphi is a potential difference rather than a field.
    const double work_base = 2.0 * b_left * s_i + gradient * s_i * s_i;
    auto work_to = [&](const double s) {
        return 2.0 * b_left * s + gradient * s * s - work_base;
    };

    // Whether the particle has the energy to reach the edge at s_b without v^2 passing
    // through zero on the way.  v^2 is quadratic in s, so when it is convex its minimum can
    // sit inside the interval even with both endpoints positive.
    auto can_reach = [&](const double s_b, double& v_sq_end) {
        v_sq_end = v_i_sqr + work_to(s_b);
        if (v_sq_end <= 0.0) { return false; }
        if (gradient > 0.0) {
            const double s_turn = -b_left / gradient;
            if ((s_turn - s_i) * (s_b - s_i) > 0.0 && (s_turn - s_b) * (s_i - s_b) > 0.0) {
                if (v_i_sqr + work_to(s_turn) <= 0.0) { return false; }
            }
        }
        return true;
    };

    // Given a reachable edge, the crossing follows from the energy relation with the
    // acceleration taken at the midpoint of the traverse, so no iteration is needed.
    // Returns false when the crossing would take longer than the sub-step allows.
    auto try_cross = [&](const double s_b, const int exit_sign, const double v_sq_end,
                         double& t_out, double& v_out) {
        v_out = exit_sign * std::sqrt(v_sq_end);
        const double mid = 0.5 * (s_i + s_b);
        const double accel_mid = (b_left * (1.0 - mid) + b_right * mid) * inv_dx;
        // pick whichever form keeps the most significant digits
        if (std::abs(v_out - v_i) > std::abs(v_out + v_i)) {
            if (accel_mid == 0.0) { return false; }
            t_out = (v_out - v_i) / accel_mid;
        } else {
            const double denom = inv_dx * (v_i + v_out);
            if (denom == 0.0) { return false; }
            t_out = 2.0 * (s_b - s_i) / denom;
        }
        return (t_out > 0.0 && t_out <= del_tau);
    };

    const int v_sign = (v_i > 0) - (v_i < 0);
    if (v_sign != 0) {
        double v_sq_end, t_cross, v_cross;
        const double s_along = (v_sign > 0) ? 1.0 : 0.0;
        if (can_reach(s_along, v_sq_end)) {
            // Reaches the edge ahead if there is time; if not it stays in the cell and the
            // iteration below resolves where.
            if (try_cross(s_along, v_sign, v_sq_end, t_cross, v_cross)) {
                s_f = s_along; v_f = v_cross; del_tau = t_cross;
                local_boundary = int(s_along);
                return true;
            }
        } else {
            // It turns around inside the cell, so the only edge it can leave by is the one
            // behind it.  This also covers a particle starting on an edge and returning to it.
            const double s_opposite = (v_sign > 0) ? 0.0 : 1.0;
            if (can_reach(s_opposite, v_sq_end)
                && try_cross(s_opposite, -v_sign, v_sq_end, t_cross, v_cross)) {
                s_f = s_opposite; v_f = v_cross; del_tau = t_cross;
                local_boundary = int(s_opposite);
                return true;
            }
        }
    }

    // Staying inside the cell: solve for the end position by Picard on the midpoint.
    double s_f_prev = s_i;
    double d_half = s_i;
    double q_dphi = b_left * (1.0 - d_half) + b_right * d_half;
    double accel = q_dphi * inv_dx;
    v_f = v_i + accel * del_tau;
    s_f = s_i + 0.5 * (v_i + v_f) * del_tau * inv_dx;
    int iteration = 0;
    while (std::abs(s_f - s_f_prev) > picard_tolerance) {
        if (++iteration > max_picard_iterations) { break; }
        s_f_prev = s_f;
        d_half = 0.5 * (s_i + s_f_prev);
        q_dphi = b_left * (1.0 - d_half) + b_right * d_half;
        accel = q_dphi * inv_dx;
        v_f = v_i + accel * del_tau;
        s_f = s_i + 0.5 * (v_i + v_f) * del_tau * inv_dx;
        if (s_f > 1.0) { s_f = 1.0; } else if (s_f < 0.0) { s_f = 0.0; }
    }
    return false;
}

} // namespace


void charged_particle::deposit_particles_quadratic(const int thread_id, std::vector<double>& work_space,
    const int left_boundary, const int right_boundary, const int number_cells) const {
    const size_t last_idx = this->number_particles[thread_id][0];
    const std::vector<double>& xi_local = this->xi[thread_id];
    for (size_t part_indx = 0; part_indx < last_idx; part_indx++){
        const double xi_p = xi_local[part_indx];
        int cell = int(xi_p);
        if (cell >= number_cells) { cell = number_cells - 1; } // particle sitting on the periodic seam
        deposit_quadratic(work_space, xi_p, cell, left_boundary, right_boundary, number_cells);
    }
}


// Push every particle over del_t with the given (time centred) field and deposit the resulting charge onto the cell centred nodes.  Particle state is left untouched, this is the residual evaluation driving the non-linear field solve.
void charged_particle::ES_push_deposit_ICIC(const int thread_id, double del_t, std::vector<double>& work_space,
    const std::vector<CIC_cell_coefficients>& cells,
    const int left_boundary, const int right_boundary, const int number_cells) {
    const std::vector<double>& xi_local = this->xi[thread_id];
    const std::vector<double>& v_x_local = this->v_x[thread_id];
    const double t_tol = del_t * 1e-10;

    auto push_range = [&](const size_t start_indx, const size_t end_indx, const std::vector<double>* del_t_array) {
        for (size_t part_indx = start_indx; part_indx < end_indx; part_indx++){
            double xi_i = xi_local[part_indx];
            double v_x_i = v_x_local[part_indx];
            const double del_t_local = (del_t_array == nullptr) ? del_t : (*del_t_array)[part_indx - start_indx];
            // nudge off an exact node so the cell is picked along the direction of travel
            int cell = int(xi_i + ((v_x_i > 0) - (v_x_i < 0)) * 1e-12);
            if (cell < 0) { cell = 0; }
            if (cell > number_cells-1) { cell = number_cells-1; }
            double del_tau = del_t_local;
            double time_passed = 0.0;
            double xi_f = xi_i, v_x_f = v_x_i;
            bool del_part = false;
            int sub_step_count = 0;
            while (del_tau > t_tol) {
                if (++sub_step_count > max_sub_steps) { break; }
                if (del_tau > cells[cell].del_tau_min) { del_tau = cells[cell].del_tau_min; }
                double s_f;
                int local_boundary = 0;
                const bool future_boundary_bool = CIC_sub_step(xi_i - double(cell), v_x_i, del_tau, s_f, v_x_f,
                    cells[cell], local_boundary);
                xi_f = double(cell) + s_f;
                if (future_boundary_bool) {
                    const int xi_boundary = cell + local_boundary;
                    if (xi_boundary == 0) {
                        switch (left_boundary){
                            case 1:
                            case 4:
                                del_part = true;
                                break;
                            case 2:
                                if (v_x_f < 0.0) { v_x_f = -v_x_f; }
                                break;
                            case 3:
                                xi_f = double(number_cells);
                                cell = number_cells-1;
                                break;
                        }
                    } else if (xi_boundary == number_cells) {
                        switch (right_boundary){
                            case 1:
                            case 4:
                                del_part = true;
                                break;
                            case 2:
                                if (v_x_f > 0.0) { v_x_f = -v_x_f; }
                                break;
                            case 3:
                                xi_f = 0.0;
                                cell = 0;
                                break;
                        }
                    } else {
                        cell = cell + ((v_x_f > 0) - (v_x_f < 0));
                    }
                    if (del_part) {
                        break;
                    }
                }
                time_passed += del_tau;
                del_tau = del_t_local - time_passed;
                xi_i = xi_f;
                v_x_i = v_x_f;
            }
            if (!del_part) {
                deposit_quadratic(work_space, xi_f, cell, left_boundary, right_boundary, number_cells);
            }
        }
    };

    size_t last_idx = this->number_particles[thread_id][0];
    push_range(0, last_idx, nullptr);
    for (int inj_indx = 0; inj_indx < this->number_unique_injections; inj_indx++) {
        const size_t start_indx = last_idx;
        last_idx = start_indx + this->number_particles_injected[thread_id][inj_indx];
        push_range(start_indx, last_idx, &this->time_step_injected[thread_id][inj_indx]);
    }
}


// Final push over del_t once the field has converged.  Updates particle state, removes particles absorbed at the walls and accumulates the wall loss diagnostics.
void charged_particle::ES_push_ICIC(const int thread_id, double del_t, int& number_sub_steps,
    const std::vector<CIC_cell_coefficients>& cells,
    const int left_boundary, const int right_boundary, const int number_cells) {
    std::vector<double>& xi_local = this->xi[thread_id];
    std::vector<double>& v_x_local = this->v_x[thread_id];
    const bool use_vy = (this->number_velocity_coordinates > 1);
    const bool use_vz = (this->number_velocity_coordinates > 2);
    const bool use_y = (this->number_space_coordinates > 1);
    const bool use_z = (this->number_space_coordinates > 2);
    const double t_tol = del_t * 1e-10;
    size_t space_delete = 0;

    auto push_range = [&](const size_t start_indx, const size_t end_indx, const std::vector<double>* del_t_array) {
        for (size_t part_indx = start_indx; part_indx < end_indx; part_indx++){
            double xi_i = xi_local[part_indx];
            double v_x_i = v_x_local[part_indx];
            const double v_y = use_vy ? this->v_y[thread_id][part_indx] : 0.0;
            const double v_z = use_vz ? this->v_z[thread_id][part_indx] : 0.0;
            const double del_t_local = (del_t_array == nullptr) ? del_t : (*del_t_array)[part_indx - start_indx];
            // nudge off an exact node so the cell is picked along the direction of travel
            int cell = int(xi_i + ((v_x_i > 0) - (v_x_i < 0)) * 1e-12);
            if (cell < 0) { cell = 0; }
            if (cell > number_cells-1) { cell = number_cells-1; }
            double del_tau = del_t_local;
            double time_passed = 0.0;
            double xi_f = xi_i, v_x_f = v_x_i;
            bool del_part = false;
            int sub_step_count = 0;
            while (del_tau > t_tol) {
                if (++sub_step_count > max_sub_steps) { break; }
                number_sub_steps++;
                if (del_tau > cells[cell].del_tau_min) { del_tau = cells[cell].del_tau_min; }
                double s_f;
                int local_boundary = 0;
                const bool future_boundary_bool = CIC_sub_step(xi_i - double(cell), v_x_i, del_tau, s_f, v_x_f,
                    cells[cell], local_boundary);
                xi_f = double(cell) + s_f;
                if (future_boundary_bool) {
                    const int xi_boundary = cell + local_boundary;
                    if (xi_boundary == 0) {
                        switch (left_boundary){
                            case 1:
                            case 4:
                                this->energy_loss[thread_id][0] += v_x_f*v_x_f + v_y*v_y + v_z*v_z;
                                this->wall_loss[thread_id][0]++;
                                this->momentum_loss[thread_id][0][0] += v_x_f;
                                this->momentum_loss[thread_id][0][1] += v_y;
                                this->momentum_loss[thread_id][0][2] += v_z;
                                del_part = true;
                                break;
                            case 2:
                                if (v_x_f < 0.0) { v_x_f = -v_x_f; }
                                break;
                            case 3:
                                xi_f = double(number_cells);
                                cell = number_cells-1;
                                break;
                        }
                    } else if (xi_boundary == number_cells) {
                        switch (right_boundary){
                            case 1:
                            case 4:
                                this->energy_loss[thread_id][1] += v_x_f*v_x_f + v_y*v_y + v_z*v_z;
                                this->wall_loss[thread_id][1]++;
                                this->momentum_loss[thread_id][1][0] += v_x_f;
                                this->momentum_loss[thread_id][1][1] += v_y;
                                this->momentum_loss[thread_id][1][2] += v_z;
                                del_part = true;
                                break;
                            case 2:
                                if (v_x_f > 0.0) { v_x_f = -v_x_f; }
                                break;
                            case 3:
                                xi_f = 0.0;
                                cell = 0;
                                break;
                        }
                    } else {
                        cell = cell + ((v_x_f > 0) - (v_x_f < 0));
                    }
                    if (del_part) {
                        break;
                    }
                }
                time_passed += del_tau;
                del_tau = del_t_local - time_passed;
                xi_i = xi_f;
                v_x_i = v_x_f;
            }
            if (!del_part) {
                const size_t new_idx = part_indx - space_delete;
                xi_local[new_idx] = xi_f;
                v_x_local[new_idx] = v_x_f;
                if (use_vy) { this->v_y[thread_id][new_idx] = v_y; }
                if (use_vz) { this->v_z[thread_id][new_idx] = v_z; }
                if (use_y) { this->y[thread_id][new_idx] = this->y[thread_id][part_indx]; }
                if (use_z) { this->z[thread_id][new_idx] = this->z[thread_id][part_indx]; }
            } else {
                space_delete++;
            }
        }
    };

    size_t last_idx = this->number_particles[thread_id][0];
    push_range(0, last_idx, nullptr);
    for (int inj_indx = 0; inj_indx < this->number_unique_injections; inj_indx++) {
        const size_t start_indx = last_idx;
        last_idx = start_indx + this->number_particles_injected[thread_id][inj_indx];
        push_range(start_indx, last_idx, &this->time_step_injected[thread_id][inj_indx]);
    }
    this->number_particles[thread_id][0] = (last_idx - space_delete);
    this->number_collidable_particles[thread_id][0] = (last_idx - space_delete);
}
