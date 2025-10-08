#pragma once
#include <vector>
#include <functional>
#include <string>



class non_linear_solver { // Base class for non-linear solvers
public:

    double eps_r, eps_a;
    int max_iterations;
    int number_iterations, number_unknowns;
    size_t accum_iterations_count;
    double solver_time, accum_residual_norm;
    virtual ~non_linear_solver() = default;

    virtual void initialize_diagnostic_files(const std::string& filename) const; // initialize diagnostic files
    virtual void write_diagnostics(const std::string& filename) const; // initialize diagnostic files
    virtual void solve(std::vector<double>& x_result, // x_result first with initial guess, then pass actual result
        const std::function<void(std::vector<double>&)>& fixed_point_function) = 0; 

    virtual void print_out() const = 0;
    
};

void read_non_linear_solver_inputs(const std::string& filename, std::vector<int>& int_params, std::vector<double>& double_params);