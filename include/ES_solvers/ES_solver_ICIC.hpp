
#pragma once
#include <vector>
#include "solvers/poisson_solver_1D.hpp"
#include "solvers/poisson_solver_1D_CIC_tridiag.hpp"
#include "non_linear_solvers/AA_solver.hpp"
#include "ES_solvers/ES_solver.hpp"
#include "domain/domain.hpp"

// Implicit, energy and charge conserving scheme with cloud-in-cell interpolation.
//
// Grid layout differs from the I-NGP scheme: the potential and charge nodes sit at the
// cell centers (number_cells of them) while the electric field lives on the cell edges
// (number_nodes of them).  A particle's cell is therefore bounded by field nodes, so the
// field it sees varies linearly and continuously as it moves, and the charge it deposits
// uses the quadratic B-spline which is the integral of that linear weighting.
class ES_solver_ICIC : public ES_solver {

public:
    std::vector<double> phi_past; // potential at the start of the time step
    std::unique_ptr<non_linear_solver> implicit_solver;
    bool smoothing;

    // geometry helpers held in the solver so the mover does not need the domain
    std::vector<double> dx_cells;        // cell widths, size number_cells
    std::vector<double> center_distance; // distance between neighbouring cell centers, size number_cells-1
    double half_dx_left, half_dx_right;  // wall to nearest cell center
    // Length each field node's value spans, so that E_field[k] * node_scale[k] recovers the
    // potential difference across the node.  The mover needs that difference rather than the
    // field itself: it converts to logical coordinates with the cell width, and the two length
    // factors have to cancel for the scheme to conserve energy on a non-uniform grid.
    std::vector<double> node_scale;
    double wrap_distance;                // across the periodic seam
    double left_boundary_potential, right_boundary_potential; // includes the RF drive
    // wall potentials at the start of the step, so a driven wall can be time centred
    double left_boundary_potential_past, right_boundary_potential_past;
    std::vector<std::vector<double>> particle_work_space; // per thread deposit scratch

    // per species mover coefficients, rebuilt every field evaluation
    std::vector<std::vector<double>> accel_node;   // q/m * E at each field node
    std::vector<std::vector<double>> del_tau_min;  // largest sub-step allowed in each cell

    ES_solver_ICIC(const domain& world);
    void initialize_diagnostic_files(const std::string& filename) override;
    void write_diagnostics(const std::string& dir_name, int diag_number) override;
    void print_out() override;
    void make_EField(const domain& world) override;
    void solve_potential(double current_time, const domain& world) override;
    void solve_field_energy(const domain& world) override;
    void write_particle_densities(const std::string file_path, const std::string filename, std::vector<charged_particle>& particle_list, const domain& world) const override;
    void get_diagnostics(const domain& world, std::vector<charged_particle>& particle_list) override;
    void deposit_charge_density(const domain& world, std::vector<charged_particle>& particle_list, int thread_id) override;
    void push_particles(int thread_id, double del_t, std::vector<charged_particle>& particle_list, const domain& world) override;
    void integrate_time_step(int thread_id, double del_t, double current_time, const domain& world, std::vector<charged_particle>& particle_list) override;

private:
    void set_boundary_potentials(double current_time, const domain& world);
    void build_mover_coefficients(std::vector<charged_particle>& particle_list, const domain& world);
};
