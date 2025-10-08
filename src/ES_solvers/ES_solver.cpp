
#include <vector>
#include <omp.h>
#include <cmath>
#include "ES_solvers/ES_solver.hpp"
#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include "ES_solvers/ES_solver_EC.hpp"
#include "ES_solvers/ES_solver_MC.hpp"
#include "ES_solvers/ES_solver_INGP.hpp"
#include "globals/mpi_vars.hpp"
#include "globals/constants.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include "globals/write_functions.hpp"
#include <iomanip>


void ES_solver::set_phi(double left_voltage, double right_voltage, double RF_frequency, int left_boundary, int right_boundary) {
    this->left_voltage = left_voltage;
    this->right_voltage = right_voltage;
    this->RF_rad_frequency = 0.0; // Convert to radians
    this->RF_half_amplitude = 0.0;
    if (left_boundary != 1 && left_boundary != 4) {
        this->left_voltage = 0.0; // Set left boundary voltage
    } 
    if (right_boundary != 1 && right_boundary != 4) {
        this->right_voltage = 0.0; // Set left boundary voltage to zero
    } 
    if (left_boundary == 4) {
        this->RF_half_amplitude = this->left_voltage; // Set RF half amplitude for left boundary
        this->RF_rad_frequency = 2.0 * RF_frequency * M_PI;
        this->left_voltage = 0.0;
    } else if (right_boundary == 4) {
        this->RF_half_amplitude = this->right_voltage; // Set RF half amplitude for right boundary
        this->RF_rad_frequency = 2.0 * RF_frequency * M_PI;
        this->right_voltage = 0.0;
    } 
    this->phi[0] = this->left_voltage; // Set left boundary voltage in phi vector
    this->phi[this->phi.size()-1] = this->right_voltage; // Set right boundary voltage in phi vector
    
}

void ES_solver::initialize_diagnostic_files(const std::string& filename) {
    // Open file
    if (mpi_vars::mpi_rank == 0) {
        std::ofstream file(filename + "/phi/parameters.dat");
        if (!file) {
            std::cerr << "Error opening file for domain \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Write header (optional)
        file << "RF_rad_frequency, RF_half_amplitude, left_voltage, right_voltage \n";

        file << std::scientific << std::setprecision(8);
        file << this->RF_rad_frequency << "\t"
        << this->RF_half_amplitude << "\t"
        << this->left_voltage << "\t"
        << this->right_voltage <<
        "\n";

        file.close();

        file.open(filename + "/field_diagnostics.dat");
        if (!file) {
            std::cerr << "Error opening file for domain \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Write header (optional)
        file << "field energy (J/m^2),  gauss error \n";

        file.close();
        

    }
}

void ES_solver::solve_field_energy(const domain& world) {
    double sum = 0.0;
    // double sum_other = 0.0;
    // // rho in computational space, so already have rho * delta_x, don't need take into consideration different delta_x non-uniform
    // if (world.left_boundary_condition == 2) {
    //     sum += this->phi[0] * this->rho[0]; // only half rho considered and only getting energy in half cell
    // }
    // if (world.right_boundary_condition == 2) {
    //     sum += this->phi[world.number_cells] * this->rho[world.number_cells];
    // }
    // for (int i = 1; i < world.number_cells; i++) {
    //     sum += this->phi[i] * this->rho[i];
    // }
    if (world.domain_type == 0) {
        // uniform domain
        double dx = world.min_dx; // Cell size for uniform domain
        for (int i = 0; i < world.number_cells; i++) {
            sum += (this->phi[i] - this->phi[i+1])*(this->phi[i] - this->phi[i+1])  / dx; // multiply by cell size
        }
    } else if (world.domain_type == 1) {
        // non-uniform domain
        const std::vector<double>& dx = world.dx_dxi; // Cell size for non-uniform domain
        for (int i = 0; i < world.number_cells; ++i) {
            sum += (this->phi[i] - this->phi[i+1])*(this->phi[i] - this->phi[i+1])  / dx[i]; // multiply by cell size
        }
    }
    this->total_field_energy = 0.5 * constants::epsilon_0 * sum;; // J/m^2
}

void ES_solver::get_diagnostics(const domain& world, std::vector<charged_particle>& particle_list) {
    this->solve_field_energy(world); // Calculate total field energy
    std::fill(this->rho.begin(), this->rho.end(), 0.0); // Reset charge density
    for (int part_num = 0; part_num < particle_list.size(); part_num++){
        for (int i = 0; i < world.number_nodes; i++) {
            this->rho[i] += particle_list[part_num].density[i] * particle_list[part_num].q_times_wp; // Accumulate charge density from all particles
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, this->rho.data(), world.number_nodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // Synchronize charge density across all processes
    std::vector<double> source_term(world.number_nodes, 0.0);
    if (world.left_boundary_condition == 2) {
        source_term[0] = -this->rho[0] /constants::epsilon_0; // Change boundary phi
    } else {
        source_term[0] = this->phi[0]; // Change boundary phi
    } 

    if (world.right_boundary_condition == 2) {
        source_term[world.number_cells] = -this->rho[world.number_cells] / constants::epsilon_0; // Change boundary phi
    } else {
        source_term[world.number_cells] = this->phi[world.number_cells]; // Change boundary phi
    } 

    for (int i = 1; i < world.number_cells; ++i) {
        source_term[i] = -this->rho[i] /constants::epsilon_0; // Set right-hand side of the Poisson equation
    }
    this->gauss_error = this->poisson_solver->norm_error(this->phi, source_term); // Solve the Poisson equation
}

void ES_solver::write_diagnostics(const std::string& dir_name, int diag_number) {
    if (mpi_vars::mpi_rank == 0) {
        write_vector_to_binary_file(this->phi, this->phi.size(), dir_name + "/phi/potential_" + std::to_string(diag_number) + ".dat", 0);
        std::ofstream file(dir_name + "/field_diagnostics.dat", std::ios::app);
        if (!file) {
            std::cerr << "Error opening file for domain \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Write header (optional)

        file << std::scientific << std::setprecision(8);
        file << this->total_field_energy << "\t"
        << this->gauss_error << "\n"; // Write field energy and Gauss error

        file.close();
    }
}



void ES_solver::deposit_charge_density(const domain& world, std::vector<charged_particle>& particle_list, int thread_id) {
    // Loop over all particles and deposit charge density
    int total_thread_count = omp_get_max_threads();
    int total_rho_size = world.number_nodes;
    int num_particles = particle_list.size();
    std::vector<double>& part_work_space = charged_particle::xi_sorted[thread_id];
    // local work_space to accumulate over each particle
    std::vector<double>& local_work_space = this->work_space[thread_id];
    std::fill(local_work_space.begin(), local_work_space.end(), 0.0);
    for (int i = 0; i < num_particles; ++i) {
        const charged_particle& particle = particle_list[i];
        // Set work space to 0
        std::fill(part_work_space.begin(), part_work_space.begin() + total_rho_size, 0.0);
        particle.deposit_particles_linear(thread_id, part_work_space); // Deposit charge density for each particle
        double q_time_wp = particle.q_times_wp; // Get charge density for each particle
        for (int j = 0; j < total_rho_size; ++j) {
            local_work_space[j] += part_work_space[j] * q_time_wp; // Accumulate charge density from all particles
        }
    }
    // Accumulate charge density from all threads
    #pragma omp barrier
    #pragma omp for
    for (int i = 0; i < total_rho_size; i++) {
        double sum = 0.0;
        for (int i_thread = 0; i_thread < total_thread_count; i_thread++) {
            sum += this->work_space[i_thread][i];
        }
        this->rho[i] = sum; // Set charge density for each cell
    }
    #pragma omp barrier
    #pragma omp master
    {
        MPI_Allreduce(MPI_IN_PLACE, this->rho.data(), total_rho_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // Synchronize charge density across all processes
    }
    #pragma omp barrier

}

void ES_solver::deposit_density(std::vector<charged_particle>& particle_list, int thread_id) {
    // Loop over all particles and deposit density
    
    int total_thread_count = omp_get_max_threads();
    int total_rho_size = this->rho.size();
    int num_particles = particle_list.size();
    // local work_space to accumulate over each particle
    std::vector<double>& local_work_space = this->work_space[thread_id];
    for (int i = 0; i < num_particles; ++i) {
        charged_particle& particle = particle_list[i];
        // Set work space to 0
        std::fill(local_work_space.begin(), local_work_space.end(), 0.0);
        particle.deposit_particles_linear(thread_id, local_work_space); // Deposit density for each particle
        // Accumulate density from all threads
        #pragma omp barrier
        #pragma omp for
        for (int i = 0; i < total_rho_size; i++) {
            double sum = 0.0;
            for (int i_thread = 0; i_thread < total_thread_count; i_thread++) {
                sum += this->work_space[i_thread][i];
            }
            particle.density[i] += sum; // Set charge density for each cell
        }
        #pragma omp barrier
    }

}

void ES_solver::write_particle_densities(const std::string file_path, const std::string filename, std::vector<charged_particle>& particle_list, const domain& world) const {
    // Loop over all particles and deposit density
    
    int number_nodes = world.number_nodes;
    int number_cells = world.number_cells;
    int num_particles = particle_list.size();
    for (int i = 0; i < num_particles; ++i) {
        charged_particle& particle = particle_list[i];
        // Set work space to 0
        std::vector<double>& density = particle.density;
        if (world.left_boundary_condition == 3) {
            density[0] = density[0] + density[number_cells];
            density[number_cells] = density[0];
        } else {
            // Only half volume represented
            density[0] = 2.0 * density[0];
            density[number_cells] = density[number_cells] * 2.0;
        }
        if (world.domain_type == 0) {
            // uniform
            double del_x = world.min_dx;
            for (int i = 0;i<number_nodes; i++) {
                density[i] = density[i] * particle.weight / del_x;
            }
        } else {
            // non-uniform
            density[0] = density[0] * particle.weight / world.dx_dxi[0];
            density[number_cells] = density[number_cells] * particle.weight / world.dx_dxi[number_cells-1];
            for (int i = 1;i<number_cells; i++) {
                double cell_size = 0.5 * (world.dx_dxi[i-1] + world.dx_dxi[i]);
                density[i] = density[i] * particle.weight / cell_size;
            }
        }
        if (mpi_vars::mpi_rank == 0) {write_vector_to_binary_file(density, number_nodes, file_path + "/charged_particles/" + particle.name + "/density/" + filename, 0);}

    }

}


void ES_solver::solve_potential(double current_time, const domain& world) {
    // generate the right-hand side of the Poisson equation
    double inv_epsilon_0 = 1.0 / constants::epsilon_0; // Inverse of permittivity
    int left_boundary = world.left_boundary_condition; // Get left boundary condition
    int right_boundary = world.right_boundary_condition; // Get right boundary condition
    int number_unknowns = this->poisson_solver->number_unknowns; // Number of unknowns in the system
    if (left_boundary == 2) {
        this->phi[0] = -this->rho[0] * inv_epsilon_0; // Change boundary phi
    } else if (left_boundary == 4) {
        this->phi[0] = this->RF_half_amplitude * std::sin(this->RF_rad_frequency * current_time); // Change boundary phi
    } 

    if (right_boundary == 2) {
        this->phi[number_unknowns-1] = -this->rho[number_unknowns-1] * inv_epsilon_0; // Change boundary phi
    } else if (right_boundary == 4) {
        this->phi[number_unknowns-1] = this->RF_half_amplitude * std::sin(this->RF_rad_frequency * current_time); // Change boundary phi
    } 

    for (int i = 1; i < number_unknowns-1; ++i) {
        this->phi[i] = -this->rho[i] * inv_epsilon_0; // Set right-hand side of the Poisson equation
    }

    this->poisson_solver->solve(this->phi, this->phi); // replace phi with solution
    
}

void ES_solver::make_EField(const domain& world) {
    // Calculate the electric field from the potential
    int number_nodes = world.number_nodes; // Number of cells in the domain
    int number_cells = world.number_cells;
    double inv_dx = 1.0/world.min_dx; // Cell size
    int left_boundary = world.left_boundary_condition; // Get left boundary condition
    int right_boundary = world.right_boundary_condition; // Get right boundary condition
    for (int i = 1; i < number_nodes-1; ++i) {
        this->E_field[i] = 0.5 * (this->phi[i-1] - this->phi[i+1]) * inv_dx; // Electric field calculation
    }
    if (left_boundary == 1 || left_boundary == 4) {
        // First order at boundary consistent with rho = 0
        this->E_field[0] = (this->phi[0] - this->phi[1])*inv_dx; // Electric field at left boundary
    } else if (left_boundary == 3){
        this->E_field[0] = 0.5 * (this->phi[number_cells-1] - this->phi[1]) * inv_dx; 
        this->E_field[number_cells] = this->E_field[0]; 
    }
    if (right_boundary == 1 || right_boundary == 4) {
        this->E_field[number_cells] = (this->phi[number_cells-1] - this->phi[number_cells]) * inv_dx; 
    } 
     
}



std::unique_ptr<ES_solver> read_voltage_inputs(const std::string& filename, int scheme_type, const domain& world){
    double left_voltage, right_voltage, RF_frequency;
    if (mpi_vars::mpi_rank == 0) {
        std::cout << " "  << std::endl;
        std::cout << "Reading electric potential inputs: "  << std::endl;
        std::cout << "-------------------------- "  << std::endl;
        std::string line;
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        
        std::getline(file, line);
        std::getline(file, line);
        std::getline(file, line);
        std::getline(file, line);
        std::getline(file, line);
        std::istringstream iss(line);
        iss >> left_voltage >> right_voltage;
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> RF_frequency;
        iss.clear();
        file.close();
    }
    MPI_Bcast(&left_voltage, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&right_voltage, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&RF_frequency, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    std::unique_ptr<ES_solver> es_solver;
    if (scheme_type == 0) {
        es_solver = std::make_unique<ES_solver_MC>(world);
    } else if (scheme_type == 1) {
        es_solver = std::make_unique<ES_solver_EC>(world);
    } else if (scheme_type == 2) {
        es_solver = std::make_unique<ES_solver_INGP>(world);
    } else {
        throw std::invalid_argument("Invalid scheme type for ES solver.");
    }
    
    es_solver->set_phi(left_voltage, right_voltage, RF_frequency, world.left_boundary_condition, world.right_boundary_condition);
    es_solver->print_out();
    
    return es_solver;

}



