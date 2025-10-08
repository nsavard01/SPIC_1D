
#pragma once
#include "solvers/poisson_solver_1D.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include <string>
#include <memory>

class poisson_solver_1D_tridiag : public poisson_solver_1D {
public:
    std::vector<double> diagonal, upper, lower; // diagonal and upper matrix elements
    poisson_solver_1D_tridiag(const domain& world);
    void solve(std::vector<double>& solution, std::vector<double>& source_term) override;
    double norm_error(const std::vector<double>& solution, const std::vector<double>& source_term) const override; // pure virtual function to solve the system
};

