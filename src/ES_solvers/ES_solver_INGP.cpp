
#include <vector>
#include <omp.h>
#include "ES_solvers/ES_solver_INGP.hpp"
#include "non_linear_solvers/non_linear_solver.hpp"
#include "solvers/poisson_solver_1D_tridiag.hpp"
#include "globals/mpi_vars.hpp"
#include "globals/constants.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "globals/write_functions.hpp"

ES_solver_INGP::ES_solver_INGP(const domain& world) {
    this->phi.resize(world.number_nodes, 0.0);
    this->phi_past.resize(world.number_nodes, 0.0);
    this->rho.resize(world.number_nodes, 0.0);
    this->J.resize(world.number_cells, 0.0);
    this->E_field.resize(world.number_cells, 0.0);
    int number_threads = omp_get_max_threads();
    this->work_space.resize(number_threads);
    for (int i = 0; i < number_threads; i++) {
        this->work_space[i].resize(world.number_nodes, 0.0);
    }
    this->poisson_solver = std::make_unique<poisson_solver_1D_tridiag>(world);
    std::vector<double> double_params(3);
    std::vector<int> int_params(3);
    read_non_linear_solver_inputs("../inputs/implicit_solver.inp", int_params, double_params); // Read non-linear solver inputs
    MPI_Bcast(double_params.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Broadcast double parameters
    MPI_Bcast(int_params.data(), 3, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast integer parameters
    this->implicit_solver = std::make_unique<AA_solver>(double_params[0], double_params[1], double_params[2], int_params[2], int_params[1], world.number_nodes);
    this->implicit_solver->print_out(); // Print implicit solver parameters
    int flag;
    if (mpi_vars::mpi_rank == 0) {
        std::string line;
        std::ifstream file("../inputs/geometry.inp");
        if (!file) {
            std::cerr << "Error: Unable to open file " << "../inputs/geometry.inp" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        
        std::getline(file, line);
        std::getline(file, line);
        std::getline(file, line);
        std::getline(file, line);
        std::getline(file, line);
        std::getline(file, line);
        std::getline(file, line);
        std::istringstream iss(line);
        iss.str(line);
        iss >> flag;
        iss.clear();
        file.close();
    }
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast smoothing flag
    this->smoothing = (flag == 1); // Set smoothing flag based on input
}

void ES_solver_INGP::initialize_diagnostic_files(const std::string& filename) {
    // Open file
    if (mpi_vars::mpi_rank == 0) {
        std::ofstream file(filename + "/phi/parameters.dat");
        if (!file) {
            std::cerr << "Error opening file for domain \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Write header (optional)
        file << "RF_rad_frequency, RF_half_amplitude, left_voltage, right_voltage smoothing \n";

        file << std::scientific << std::setprecision(8);
        file << this->RF_rad_frequency << "\t"
        << this->RF_half_amplitude << "\t"
        << this->left_voltage << "\t"
        << this->right_voltage << "\t"
        << (this->smoothing ? 1 : 0) << // Write smoothing flag
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
        
        this->implicit_solver->initialize_diagnostic_files(filename);
    }
}

void ES_solver_INGP::write_diagnostics(const std::string& dir_name, int diag_number) {
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
        this->implicit_solver->write_diagnostics(dir_name); // Write implicit solver diagnostics
    }
}

void ES_solver_INGP::print_out() {
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "ES_solver_INGP: " << std::endl;
        std::cout << "-------------------------- " << std::endl;
        std::cout << "Number of phi nodes: " << this->phi.size() << std::endl;
        std::cout << "Number of field nodes: " << this->E_field.size() << std::endl;
        std::cout << "Left voltage: " << this->left_voltage << std::endl;
        std::cout << "Right voltage: " << this->right_voltage << std::endl;
        std::cout << "RF frequency: " << this->RF_rad_frequency / (2.0 * M_PI) << std::endl;
        if (this->RF_half_amplitude != 0.0) {
            std::cout << "RF half amplitude " << this->RF_half_amplitude << std::endl;
        } else {
            std::cout << "No RF set." << std::endl;
        }
        std::cout << "Smoothing: " << (this->smoothing ? "Enabled" : "Disabled") << std::endl;
        std::cout << "-------------------------- " << std::endl;
    }
}

inline void smooth_charge_density(std::vector<double>& rho, const domain& world) {
    // Smooth charge density using binomial smoothing
    std::vector<double> rho_copy = rho; // Copy charge density for smoothing
    int left_boundary = world.left_boundary_condition; // Get left boundary condition
    int right_boundary = world.right_boundary_condition; // Get right boundary condition
    double rho_edge;
    // binomial smoothing
    if (left_boundary == 1 || left_boundary == 4) {
        // rho is 0 at edge
        rho_edge = 0.0;
        rho[0] = 0.0;
    } else if (left_boundary == 2){
        rho_edge = 2.0 * rho_copy[0];
        // rho_copy[0] is half volume so multiply by 2
        rho[0] = 0.125 * (2.0 * rho_edge + 2.0 * rho_copy[1]); // stencil assumes half rho at edge
    } else {
        rho_edge = rho_copy[0] + rho_copy[world.number_cells]; // rho_copy[0] is half volume so add full volume at edge
        rho[0] = 0.25 * (rho_copy[world.number_cells-1] + 2.0 * rho_edge + rho_copy[1]);
        rho[world.number_cells] = rho[0];
    }
    rho[1] = 0.25 * (rho_edge + 2.0 * rho_copy[1] + rho_copy[2]); // Smooth first cell
    for (int i = 2; i < world.number_cells - 1; i++) {
        rho[i] = 0.25 * (rho_copy[i-1] + 2.0 * rho_copy[i] + rho_copy[i+1]); // Smooth interior cells
    }
    if (right_boundary == 1 || right_boundary == 4) {
        // rho is 0 at edge
        rho_edge = 0.0;
        rho[world.number_cells] = 0.0;
    } else if (right_boundary == 2){
        rho_edge = 2.0 * rho_copy[world.number_cells];
        // rho_copy[0] is half volume so multiply by 2
        rho[world.number_cells] = 0.125 * (2.0 * rho_edge + 2.0 * rho_copy[world.number_cells-1]); // stencil assumes half rho at edge
    } // already taken care of with left_boundary == 3
    rho[world.number_cells-1] = 0.25 * (rho_edge + 2.0 * rho_copy[world.number_cells-1] + rho_copy[world.number_cells-2]); // Smooth last cell
        
}

void ES_solver_INGP::deposit_charge_density(const domain& world, std::vector<charged_particle>& particle_list, int thread_id) {
    // Loop over all particles and deposit charge density
    
    int total_thread_count = omp_get_max_threads();
    int total_rho_size = world.number_nodes; // Total number of nodes in the domain
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
        if (this->smoothing) {
            smooth_charge_density(this->rho, world); // Smooth charge density if enabled 
        }
    }
    #pragma omp barrier

}

inline void smooth_field(std::vector<double>& E_field, const domain& world) {
    // Smooth electric field using binomial smoothing
    std::vector<double> E_field_copy = E_field; // Copy electric field for smoothing
    int left_boundary = world.left_boundary_condition; // Get left boundary condition
    int right_boundary = world.right_boundary_condition; // Get right boundary condition
    // binomial smoothing
    if (left_boundary == 1 || left_boundary == 4) {
        // E_field is 0 at edge
        E_field[0] = 0.25 * (3.0 * E_field_copy[0] + E_field_copy[1]);
    } else if (left_boundary == 2){
        E_field[0] = 0.25 * (E_field_copy[0] + E_field_copy[1]);
    } else {
        E_field[0] = 0.25 * (2.0 * E_field_copy[0] + E_field_copy[1] + E_field_copy[world.number_cells-1]); // Smooth first cell
    }
    for (int i = 1; i < world.number_cells - 1; i++) {
        E_field[i] = 0.25 * (E_field_copy[i-1] + 2.0 * E_field_copy[i] + E_field_copy[i+1]);
    }
    if (right_boundary == 1 || right_boundary == 4) {
        // E_field is 0 at edge
        E_field[world.number_cells-1] = 0.25 * (3.0 * E_field_copy[world.number_cells-1] + E_field_copy[world.number_cells-2]);
    } else if (right_boundary == 2){
        E_field[world.number_cells-1] = 0.25 * (E_field_copy[world.number_cells-1] + E_field_copy[world.number_cells-2]);
    } else {
        E_field[world.number_cells-1] = 0.25 * (2.0 * E_field_copy[world.number_cells-1] + E_field_copy[0] + E_field_copy[world.number_cells-2]); // Smooth first cell
    }
}


void ES_solver_INGP::make_EField(const domain& world) {
    // Calculate the electric field from the potential
    int number_cells = world.number_cells; // Number of cells in the domain
    for (int i = 0; i < number_cells; ++i) {
        this->E_field[i] = 0.5 * (this->phi[i] + this->phi_past[i] - this->phi[i+1] - this->phi_past[i+1]); // Electric field calculation
    }
    if (this->smoothing) {
        smooth_field(this->E_field, world); // Smooth electric field if smoothing is enabled
    }
    if (world.domain_type == 0) {
        double inv_dx = 1.0/world.min_dx; // Cell size for uniform domain
        inv_dx = 1.0/world.min_dx; // Cell size for uniform domain
        for (int i = 0; i < number_cells; ++i) {
            this->E_field[i] = this->E_field[i] * inv_dx; // Electric field calculation
        }
    } else if (world.domain_type == 1) {
        const std::vector<double>& dx = world.dx_dxi; // Cell size for non-uniform domain
        for (int i = 0; i < number_cells; ++i) {
            this->E_field[i] = this->E_field[i]/ dx[i]; // Electric field calculation
        }
    }

}

void ES_solver_INGP::get_diagnostics(const domain& world, std::vector<charged_particle>& particle_list) {
    this->solve_field_energy(world); // Calculate total field energy
    std::fill(this->rho.begin(), this->rho.end(), 0.0); // Reset charge density
    for (int part_num = 0; part_num < particle_list.size(); part_num++){
        for (int i = 0; i < world.number_nodes; i++) {
            this->rho[i] += particle_list[part_num].density[i] * particle_list[part_num].q_times_wp; // Accumulate charge density from all particles
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, this->rho.data(), world.number_nodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // Synchronize charge density across all processes
    if (this->smoothing) {
        smooth_charge_density(this->rho, world); // Smooth charge density if enabled 
    }
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



void ES_solver_INGP::push_particles(const int thread_id, double del_t, std::vector<charged_particle>& particle_list, const domain& world){
    // Loop over all particles and push them to the grid
    int num_particles = particle_list.size();
    int left_boundary = world.left_boundary_condition; // Get left boundary condition
    int right_boundary = world.right_boundary_condition; // Get right boundary condition
    int number_cells = world.number_cells; // Number of cells in the domain
    int domain_type = world.domain_type;
    
    std::vector<double>& local_work_space = this->work_space[thread_id];
    std::vector<double>& part_work_space = charged_particle::xi_sorted[thread_id];
    std::fill(local_work_space.begin(), local_work_space.end(), 0.0); // Reset work space for this thread
    if (domain_type == 0) {
        double inv_dx = 1.0 / world.min_dx; // Cell size
        for (int i = 0; i < num_particles; ++i) {
            charged_particle& particle = particle_list[i];
            std::fill(part_work_space.begin(), part_work_space.begin() + world.number_nodes, 0.0); // Reset work space for this particle
            particle.ES_push_deposit_INGP_uniform(thread_id, del_t, this->E_field, part_work_space, 
                inv_dx, left_boundary, right_boundary, number_cells); // Push particles to the grid
            for (int i = 0; i < world.number_nodes; i++) {
                local_work_space[i] += part_work_space[i] * particle.q_times_wp;
            }
        }
    } else {
        for (int i = 0; i < num_particles; ++i) {
            charged_particle& particle = particle_list[i];
            std::fill(part_work_space.begin(), part_work_space.begin() + world.number_nodes, 0.0); // Reset work space for this particle
            particle.ES_push_deposit_INGP_non_uniform(thread_id, del_t, this->E_field, part_work_space, 
                world.dx_dxi, left_boundary, right_boundary, number_cells); // Push particles to the grid
            for (int i = 0; i < world.number_nodes; i++) {
                local_work_space[i] += part_work_space[i] * particle.q_times_wp;
            }
        }
    }
}

void ES_solver_INGP::write_particle_densities(const std::string file_path, const std::string filename, std::vector<charged_particle>& particle_list, const domain& world) const {
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
        // MPI_Allreduce(MPI_IN_PLACE, density.data(), number_nodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (this->smoothing) {
            std::vector<double> density_copy = density; // Copy density for smoothing
            if (world.left_boundary_condition == 3) {
                // density is 0 at edge
                density[0] = 0.25 * (2.0 * density_copy[0] + density_copy[1] + density_copy[number_cells]);
            } else {
                // Even for dirichlet assume symmetry at edge
                density[0] = 0.25 * (2.0 * density_copy[0] + 2.0 * density_copy[1]); // Smooth first cell
            }
            for (int j = 1; j < number_cells; j++) {
                density[j] = 0.25 * (density_copy[j-1] + 2.0 * density_copy[j] + density_copy[j+1]); // Smooth interior cells
            }
            if (world.right_boundary_condition == 3) {
                // density is 0 at edge
                density[number_cells] = density[0];
            } else {
                // Even for dirichlet assume symmetry at edge
                density[number_cells] = 0.25 * (2.0 * density_copy[number_cells] + 2.0 * density_copy[number_cells]); // Smooth first cell
            }
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

void ES_solver_INGP::integrate_time_step(const int thread_id, double del_t, double current_time, const domain& world, std::vector<charged_particle>& particle_list) {

    auto integral_function = [&](std::vector<double>& res_output) {
        int total_thread_count = omp_get_max_threads();
        double part_timer_start;
        #pragma omp master
        {   
            this->make_EField(world);
            part_timer_start = MPI_Wtime();
        }
        #pragma omp barrier
        this->push_particles(thread_id, del_t, particle_list, world);
        #pragma omp barrier
        #pragma omp for
        for (int i = 0; i < world.number_nodes; i++) {
            double sum = 0.0;
            for (int i_thread = 0; i_thread < total_thread_count; i_thread++) {
                sum += this->work_space[i_thread][i];
            }
            this->rho[i] = sum; // Set charge density for each cell
        }
        #pragma omp barrier
        #pragma omp master
        {
            MPI_Allreduce(MPI_IN_PLACE, this->rho.data(), world.number_nodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // Synchronize charge density across all processes
            if (this->smoothing) {
                smooth_charge_density(this->rho, world); // Smooth charge density if enabled 
            }
        }
        #pragma omp barrier
        #pragma omp master
        {
            double end_time = MPI_Wtime();
            this->particle_timer += (end_time - part_timer_start);
            double start_time = MPI_Wtime();
            this->solve_potential(current_time + del_t, world);
            end_time = MPI_Wtime();
            this->potential_timer += (end_time - start_time);
        }
    };
    // double KE_initial = 0.0, PE_initial=0.0, KE_final=0.0, PE_final=0.0;
    #pragma omp master
    {
        this->phi_past = this->phi; // Copy current potential to future potential
        this->particle_timer = 0.0;
        this->potential_timer = 0.0;
        // PE_initial = this->total_field_energy; // Store initial potential energy
        // KE_initial = 0.0; // Initialize initial kinetic energy
        // for (const auto& particle : particle_list) {
        //     double sum = 0.0;
        //     for (int i = 0; i < particle.number_velocity_coordinates; ++i) {
        //         sum += particle.total_sum_v_square[i]; // Accumulate velocity
        //     }
        //     KE_initial += sum * 0.5 * particle.weight * particle.mass; // Calculate initial kinetic energy
        // }
        // if (mpi_vars::mpi_rank == 0) {
        //     std::cout << "ES_solver_INGP: Integrating time step at time " << current_time << " with del_t = " << del_t << std::endl;
        //     std::cout << "Initial KE = " << KE_initial << ", Initial PE = " << PE_initial << std::endl;
        //     std::cout << "Initial total Energy: " << KE_initial + PE_initial << std::endl;
        // }
    }
    #pragma omp barrier
    this->implicit_solver->solve(this->phi, integral_function); // Solve the non-linear system
    #pragma omp master
    {
        this->make_EField(world); // Calculate the electric field from the potential
    }
    #pragma omp barrier
    int num_particles = particle_list.size();
    int left_boundary = world.left_boundary_condition; // Get left boundary condition
    int right_boundary = world.right_boundary_condition; // Get right boundary condition
    int number_cells = world.number_cells; // Number of cells in the domain
    int domain_type = world.domain_type;
    std::vector<size_t> number_initial_particles(num_particles, 0); // Initialize number of initial particles for each particle
    std::vector<int> number_sub_steps(num_particles, 0); // Initialize number of sub-steps for each particle
    // Final push of particles
    if (domain_type == 0) {
        double inv_dx = 1.0 / world.min_dx; // Cell size
        for (int i = 0; i < num_particles; ++i) {
            charged_particle& particle = particle_list[i];
            number_initial_particles[i] = particle.number_particles[thread_id][0]; // Store initial number of particles
            particle.ES_push_INGP_uniform(thread_id, del_t, this->E_field, number_sub_steps[i], 
                inv_dx, left_boundary, right_boundary, number_cells); // Push particles to the grid
        }
    } else {
        for (int i = 0; i < num_particles; ++i) {
            charged_particle& particle = particle_list[i];
            particle.ES_push_INGP_non_uniform(thread_id, del_t, this->E_field, number_sub_steps[i], 
                world.dx_dxi, left_boundary, right_boundary, number_cells); // Push particles to the grid
        }
    }
    #pragma omp barrier
   
    
    // this->deposit_charge_density(world, particle_list, thread_id);
    // #pragma omp barrier
    // for (int part_num = 0; part_num < particle_list.size(); part_num++){
    //     particle_list[part_num].get_particle_diagnostics(thread_id, world.number_cells); // Gather particle diagnostics
    // }
    // #pragma omp barrier
    // #pragma omp master
    // {
    //     this->solve_field_energy(world);
    //     PE_final = this->total_field_energy; // Store final potential energy
    //     for (int part_num = 0; part_num < particle_list.size(); part_num++){
    //         particle_list[part_num].gather_mpi(); // Gather particle diagnostics
    //     }
    //     KE_final = 0.0; // Initialize final kinetic energy
    //     for (const auto& particle : particle_list) {
    //         double sum = 0.0;
    //         for (int i = 0; i < particle.number_velocity_coordinates; ++i) {
    //             sum += particle.total_sum_v_square[i]; // Accumulate velocity
    //         }
    //         if (mpi_vars::mpi_rank == 0) {
    //             std::cout << "Number of " << particle.name << " is " << particle.total_number_particles << std::endl; // Write particle diagnostics
    //         }
    //         KE_final += 0.5 * particle.weight * particle.mass * (sum + particle.accum_wall_energy_loss[0] + particle.accum_wall_energy_loss[1]); // Calculate final kinetic energy
    //     }
    //     std::vector<double> source_term(world.number_nodes, 0.0);
    //     if (world.left_boundary_condition == 2) {
    //         source_term[0] = -this->rho[0] /constants::epsilon_0; // Change boundary phi
    //     } else {
    //         source_term[0] = this->phi[0]; // Change boundary phi
    //     } 
    
    //     if (world.right_boundary_condition == 2) {
    //         source_term[world.number_cells] = -this->rho[world.number_cells] / constants::epsilon_0; // Change boundary phi
    //     } else {
    //         source_term[world.number_cells] = this->phi[world.number_cells]; // Change boundary phi
    //     } 
    
    //     for (int i = 1; i < world.number_cells; ++i) {
    //         source_term[i] = -this->rho[i] /constants::epsilon_0; // Set right-hand side of the Poisson equation
    //     }
    //     double error = this->poisson_solver->norm_error(this->phi, source_term); // Solve the Poisson equation
    //     if (error > 1e-4) {
    //         if (mpi_vars::mpi_rank == 0) {
    //             std::cerr << "ES_solver_INGP: Norm error in Poisson solver is too high: " << error << std::endl;
    //             std::cerr << "This may indicate a problem with the solver or the input parameters." << std::endl;
    //         }
    //     }
    //     if (mpi_vars::mpi_rank == 0) {
    //         std::cout << "PE_initial = " << PE_initial << ", PE_final = " << PE_final << std::endl;
    //         std::cout << "Final KE = " << KE_final << ", Final PE = " << PE_final << std::endl;
    //         std::cout << "ES_solver_INGP: Final total Energy: " << KE_final + PE_final << std::endl;
    //         std::cout << "Difference in total energy: " << ((KE_final + PE_final) - (KE_initial + PE_initial))/(KE_initial + PE_initial) << std::endl;
    //         std::cout << "took " << this->implicit_solver->number_iterations << " iterations to converge." << std::endl;
    //         std::cout << "ES_solver_INGP: Norm error in Poisson solver: " << error << std::endl;
    //     }
    // }
    
    
    

}





