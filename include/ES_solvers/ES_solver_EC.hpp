
#pragma once
#include <vector>
#include "solvers/poisson_solver_1D.hpp"
#include "ES_solvers/ES_solver.hpp"
#include "domain/domain.hpp"

class ES_solver_EC : public ES_solver {


public:
    ES_solver_EC(const domain& world);
    void print_out() override;
    void make_EField(const domain& world) override;
    void push_particles(int thread_id, double del_t, std::vector<charged_particle>& particle_list, const domain& world) override;
    void integrate_time_step(int thread_id, double del_t, double current_time, const domain& world, std::vector<charged_particle>& particle_list) override;
    // void integrate_time(double current_time, double del_t) override;
    // void interpolate_particles_to_grid(std::vector<charged_particle>& particle_list) override;
};

