#include "solvers/poisson_solver_1D_CIC_tridiag.hpp"
#include "globals/mpi_vars.hpp"
#include <cmath>
#include <stdexcept>

// Row i integrates Poisson's equation over cell i, whose edges are the field (half) nodes:
//   (phi_{i+1} - phi_i)/h_i - (phi_i - phi_{i-1})/h_{i-1} = -rho_i / epsilon_0
// with h_i the distance between neighbouring cell centers.  At a Dirichlet wall the
// outward gradient is taken over the half cell separating the wall from the first center.
poisson_solver_1D_CIC_tridiag::poisson_solver_1D_CIC_tridiag(const domain& world) {
    this->number_unknowns = world.number_cells;
    this->diagonal.resize(this->number_unknowns, 0.0);
    this->upper.resize(this->number_unknowns-1, 0.0);
    this->lower.resize(this->number_unknowns-1, 0.0);
    this->work_space.resize(this->number_unknowns, 0.0);

    const int number_cells = world.number_cells;
    const std::vector<double>& centers = world.cell_centers;
    const std::vector<double>& nodes = world.grid_nodes;
    // distance from each wall to the nearest cell center
    const double half_dx_left = centers[0] - nodes[0];
    const double half_dx_right = nodes[number_cells] - centers[number_cells-1];

    // interior couplings, distance between neighbouring cell centers
    for (int i = 0; i < number_cells - 1; ++i) {
        double h = centers[i+1] - centers[i];
        this->upper[i] = 1.0 / h;
        this->lower[i] = 1.0 / h;
    }
    for (int i = 1; i < number_cells - 1; ++i) {
        this->diagonal[i] = -(this->lower[i-1] + this->upper[i]);
    }

    // Boundary rows.  Periodic is closed the same way the nodal (I-NGP) solver closes it,
    // by pinning the potential on the wall itself, which is a valid gauge choice for a
    // periodic domain but does assume the wall stays at the given potential.
    switch (world.left_boundary_condition) {
        case 1:
        case 4:
        case 3:
            this->diagonal[0] = -(1.0 / half_dx_left + this->upper[0]);
            break;
        case 2:
            // reflective/symmetric, no field through the wall
            this->diagonal[0] = -this->upper[0];
            break;
        default:
            throw std::invalid_argument("Invalid left boundary condition.");
    }

    switch (world.right_boundary_condition) {
        case 1:
        case 4:
        case 3:
            this->diagonal[number_cells-1] = -(1.0 / half_dx_right + this->lower[number_cells-2]);
            break;
        case 2:
            this->diagonal[number_cells-1] = -this->lower[number_cells-2];
            break;
        default:
            throw std::invalid_argument("Invalid right boundary condition.");
    }
}

void poisson_solver_1D_CIC_tridiag::solve(std::vector<double>& solution, std::vector<double>& source_term) {
    // Thomas algorithm, solution and source_term may alias
    double m;
    this->work_space[0] = this->upper[0] / this->diagonal[0];
    solution[0] = source_term[0] / this->diagonal[0];
    for (int i = 1; i < this->number_unknowns-1; i++) {
        m = this->diagonal[i] - this->lower[i-1] * this->work_space[i-1];
        this->work_space[i] = this->upper[i] / m;
        solution[i] = (source_term[i] - this->lower[i-1] * solution[i-1]) / m;
    }
    m = this->diagonal[this->number_unknowns-1] - this->lower[this->number_unknowns-2] * this->work_space[this->number_unknowns-2];
    solution[this->number_unknowns-1] = (source_term[this->number_unknowns-1] - this->lower[this->number_unknowns-2] * solution[this->number_unknowns-2]) / m;
    for (int i = this->number_unknowns-2; i >= 0; i--) {
        solution[i] = solution[i] - this->work_space[i] * solution[i+1];
    }
}

double poisson_solver_1D_CIC_tridiag::norm_error(const std::vector<double>& solution, const std::vector<double>& source_term) const {
    double error = 0.0;
    double res;
    if (source_term[0] != 0.0) {
        res = (this->diagonal[0] * solution[0] + this->upper[0] * solution[1]) / source_term[0] - 1.0;
        error += res * res;
    }
    for (int i = 1; i < this->number_unknowns-1; i++) {
        if (source_term[i] != 0.0) {
            res = (this->lower[i-1] * solution[i-1] + this->diagonal[i] * solution[i] + this->upper[i] * solution[i+1]) / source_term[i] - 1.0;
            error += res * res;
        }
    }
    if (source_term[this->number_unknowns-1] != 0.0) {
        res = (this->lower[this->number_unknowns-2] * solution[this->number_unknowns-2] + this->diagonal[this->number_unknowns-1] * solution[this->number_unknowns-1]) / source_term[this->number_unknowns-1] - 1.0;
        error += res * res;
    }
    return std::sqrt(error / static_cast<double>(this->number_unknowns));
}
