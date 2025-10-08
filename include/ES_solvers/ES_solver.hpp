
#pragma once
#include <vector>
#include <cmath>
#include "solvers/poisson_solver_1D.hpp"
#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include "globals/mpi_vars.hpp"
#include "particles/charged_particle.hpp"
#include <cmath>

class ES_solver {

public:
    std::vector<double> phi, rho, E_field; // diagonal and upper matrix elements
    double gauss_error = 0.0;
    double RF_half_amplitude, RF_rad_frequency, left_voltage, right_voltage; // RF amplitude and frequency
    std::unique_ptr<poisson_solver_1D> poisson_solver; // pointer to the Poisson solver
    std::vector<std::vector<double>> work_space;

    // diagnostics
    double total_field_energy;
    double potential_timer, particle_timer;
    
    virtual ~ES_solver() = default;
    
    std::vector<double>& get_phi() {
        return this->phi;
    };

    std::vector<double>& get_rho() {
        return this->rho;
    };
    virtual void initialize_diagnostic_files(const std::string& filename);
    virtual void write_diagnostics(const std::string& dir_name, int diag_number);

    void set_phi(double left_voltage, double right_voltage, double RF_frequency, int left_boundary, int right_boundary);

    virtual void print_out() = 0;
    virtual void write_particle_densities(const std::string file_path, const std::string filename, std::vector<charged_particle>& particle_list, const domain& world) const; // since density determinined by potential solver type
    virtual void deposit_charge_density(const domain& world, std::vector<charged_particle>& particle_list, int thread_id);
    virtual void deposit_density(std::vector<charged_particle>& particle_list, int thread_id);
    virtual void solve_potential(double current_time, const domain& world);
    virtual void solve_field_energy(const domain& world);
    virtual void make_EField(const domain& world);
    // general integration through time step which solves for fields after time step given initial fields
    // general enough that it can include non-linear processes as well (implicit)
    virtual void get_diagnostics(const domain& world, std::vector<charged_particle>& particle_list);
    virtual void integrate_time_step(const int thread_id, double del_t, double current_time, const domain& world, std::vector<charged_particle>& particle_list) = 0; 
    virtual void push_particles(const int thread_id, double del_t, std::vector<charged_particle>& particle_list, const domain& world) = 0;

};

std::unique_ptr<ES_solver> read_voltage_inputs(const std::string& filename, int scheme_type, const domain& world);

