
#pragma once
#include <vector>
#include "solvers/poisson_solver_1D.hpp"
#include "non_linear_solvers/AA_solver.hpp"
#include "ES_solvers/ES_solver.hpp"
#include "domain/domain.hpp"

class ES_solver_INGP : public ES_solver {


public:
    std::vector<double> phi_past, J; // future potential for implicit solver
    std::unique_ptr<non_linear_solver> implicit_solver; // pointer to the implicit non-linear solver
    bool smoothing;
    ES_solver_INGP(const domain& world);
    void initialize_diagnostic_files(const std::string& filename) override;
    void write_diagnostics(const std::string& dir_name, int diag_number) override;
    void print_out() override;
    void make_EField(const domain& world) override;
    void write_particle_densities(const std::string file_path, const std::string filename, std::vector<charged_particle>& particle_list, const domain& world) const override;
    virtual void get_diagnostics(const domain& world, std::vector<charged_particle>& particle_list) override; 
    void deposit_charge_density(const domain& world, std::vector<charged_particle>& particle_list, int thread_id) override;
    void push_particles(int thread_id, double del_t, std::vector<charged_particle>& particle_list, const domain& world) override;
    void integrate_time_step(int thread_id, double del_t, double current_time, const domain& world, std::vector<charged_particle>& particle_list) override;
    // void integrate_time(double current_time, double del_t) override;
    // void interpolate_particles_to_grid(std::vector<charged_particle>& particle_list) override;
};

