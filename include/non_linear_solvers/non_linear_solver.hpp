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
    // convergence bookkeeping: a step that exhausts max_iterations is accepted anyway,
    // so it leaves no other trace.  These are running values over the whole run.
    int max_iterations_used;        // most iterations any single step has needed
    size_t non_converged_count;     // steps that hit max_iterations without converging
    double worst_residual;          // largest final residual of any step
    double convergence_tolerance;   // eps_a * sqrt(N), the absolute part of the test
    int warnings_printed;
    std::vector<double> residual_history;   // residual at each iteration of the current step
    virtual ~non_linear_solver() = default;

    virtual void initialize_diagnostic_files(const std::string& filename) const; // initialize diagnostic files
    virtual void write_diagnostics(const std::string& filename) const; // initialize diagnostic files
    virtual void solve(std::vector<double>& x_result, // x_result first with initial guess, then pass actual result
        const std::function<void(std::vector<double>&)>& fixed_point_function) = 0; 

    virtual void print_out() const = 0;
    
};

void read_non_linear_solver_inputs(const std::string& filename, std::vector<int>& int_params, std::vector<double>& double_params);