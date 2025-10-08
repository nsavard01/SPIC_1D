
#pragma once
#include <vector>
#include <string>
#include <memory>

class poisson_solver_1D {

public:
    int number_unknowns; // number of unknowns in the system
    std::vector<double> work_space; // diagonal and upper matrix elements
    
    virtual ~poisson_solver_1D() = default;

    virtual void solve(std::vector<double>& solution, std::vector<double>& source_term) = 0; // pure virtual function to solve the system
    virtual double norm_error(const std::vector<double>& solution, const std::vector<double>& source_term) const = 0; // pure virtual function to solve the system
};

