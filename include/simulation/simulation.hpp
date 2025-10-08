
#pragma once

#include "globals/constants.hpp"
#include "globals/mpi_vars.hpp"
#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include "rand_gen/pcg_rng.hpp"
#include "particles/charged_particle.hpp"
#include "particles/target_particle.hpp"
#include "solvers/poisson_solver_1D_tridiag.hpp"
#include "ES_solvers/ES_solver_MC.hpp"
#include "ES_solvers/ES_solver_EC.hpp"
#include "collisions/null_collider.hpp"
#include "charged_particle_operators/charged_particle_operator.hpp"
#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <dirent.h>
#include <sys/stat.h>

class simulation {
    
public:
    // Everything needed to simulate particle in cell
    double del_t, simulation_time, simulation_start_time, averaging_time, current_time, diag_time_division, next_diag_time, last_diag_time;
    double elapsed_time, particle_time, field_time, null_collision_time, art_collision_time;
    double res_acceptance_phi, res_acceptance_density;
    int current_diag_step;
    size_t diag_step_diff, current_step; 
    double inv_plasma_freq_fraction;
    bool restarted_simulation;
    std::string save_file_folder, save_file_path;
    int number_omp_threads, scheme_type, number_diagnostics;
    std::unique_ptr<domain> world;
    std::vector<charged_particle> charged_particle_list;
    std::vector<target_particle> target_particle_list;
    std::vector<null_collider> null_collider_list;
    std::unique_ptr<ES_solver> field_solver;
    std::vector<std::unique_ptr<charged_particle_operator>> particle_operator_list;
    simulation();
    void setup();
    void initialize_diagnostic_files();
    void diagnostics(int thread_id);
    void reset_diagnostics(int thread_id);
    void run();
    void averaging();
};

// Function to remove a directory and its contents
inline void removeDirectoryContents(const std::string& dirName) {
    std::string command = "rm -r " + dirName + "/*";
    int status = system(command.c_str());
    if (status != 0) {
        std::cerr << "Error removing directory contents: " << dirName << std::endl;
        exit(1);  // Exiting since the function isn't able to clean up properly
    }
}

// Function to create a directory
inline bool createDirectory(const std::string& dirName) {
    if (mkdir(dirName.c_str(), 0777) == -1) {
        std::cerr << "Error creating directory: " << dirName << std::endl;
        return false;
    }
    return true;
}

inline bool directoryExists(const std::string& dirName) {
    DIR* dir = opendir(dirName.c_str());
    if (dir) {
        closedir(dir);
        return true;
    }
    return false;
}




