
#pragma once
#include <vector>
#include <cmath>
#include "domain/domain.hpp"
#include "globals/mpi_vars.hpp"
#include "particles/charged_particle.hpp"
#include <cmath>
#include "globals/constants.hpp"
#include <omp.h>
#include <fstream>
#include <sstream>
#include <iomanip>

class charged_particle_operator {

// general charged particle operator base
// input charged particle list, thread_id, and domain, and apply with run function using internal variables
public:

    std::vector<int> particle_indx;
    virtual ~charged_particle_operator() = default;
    virtual void print_out() = 0;
    virtual void setup_diagnostics(const std::string& dir_name, std::vector<charged_particle>& particle_list) = 0;
    virtual void write_diagnostics(const std::string& dir_name, const std::vector<charged_particle>& particle_list) = 0;
    virtual void reset_diagnostics() = 0;
    virtual void write_average_diagnostics(const std::string& dir_name, const std::vector<charged_particle>& particle_list) = 0;
    virtual void run(const int thread_id, const double current_time, const double del_t, std::vector<charged_particle>& particle_list, const domain& world) = 0;
};

std::vector<std::unique_ptr<charged_particle_operator>> read_particle_operators(const std::string& directory_path, std::vector<charged_particle>& particle_list, const domain& world);
