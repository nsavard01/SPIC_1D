
#pragma once
#include "solvers/poisson_solver_1D.hpp"
#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include <string>
#include <memory>

// Tridiagonal Poisson solver for the CIC (cell-centered) grid arrangement.
// Unlike the nodal solver, the unknowns phi_i sit at the cell centers, so there is
// no unknown sitting on the wall itself.  Boundary values are applied half a cell
// outside of the first/last unknown, which is why the boundary rows are not simply
// identity rows like they are in poisson_solver_1D_tridiag.
class poisson_solver_1D_CIC_tridiag : public poisson_solver_1D {
public:
    std::vector<double> diagonal, upper, lower; // diagonal, upper and lower matrix elements
    poisson_solver_1D_CIC_tridiag(const domain& world);
    void solve(std::vector<double>& solution, std::vector<double>& source_term) override;
    double norm_error(const std::vector<double>& solution, const std::vector<double>& source_term) const override;
};
