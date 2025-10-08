
#include <vector>
#include <omp.h>
#include "ES_solvers/ES_solver_EC.hpp"
#include "solvers/poisson_solver_1D_tridiag.hpp"
#include "globals/mpi_vars.hpp"
#include "globals/constants.hpp"

ES_solver_EC::ES_solver_EC(const domain& world) {
    this->phi.resize(world.number_nodes, 0.0);
    this->rho.resize(world.number_nodes, 0.0);
    this->E_field.resize(world.number_cells, 0.0);
    int number_threads = omp_get_max_threads();
    this->work_space.resize(number_threads);
    for (int i = 0; i < number_threads; i++) {
        this->work_space[i].resize(world.number_nodes, 0.0);
    }
    this->poisson_solver = std::make_unique<poisson_solver_1D_tridiag>(world);
}

void ES_solver_EC::print_out() {
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "ES_solver_EC: " << std::endl;
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
        std::cout << "-------------------------- " << std::endl;
    }
}


void ES_solver_EC::make_EField(const domain& world) {
    // Calculate the electric field from the potential
    int number_cells = world.number_cells; // Number of cells in the domain
    if (world.domain_type == 0) {
        double inv_dx = 1.0/world.min_dx; // Cell size for uniform domain
        inv_dx = 1.0/world.min_dx; // Cell size for uniform domain
        for (int i = 0; i < number_cells; ++i) {
            this->E_field[i] = (this->phi[i] - this->phi[i+1]) * inv_dx; // Electric field calculation
        }
    } else if (world.domain_type == 1) {
        const std::vector<double>& dx = world.dx_dxi; // Cell size for non-uniform domain
        for (int i = 0; i < number_cells; ++i) {
            this->E_field[i] = (this->phi[i] - this->phi[i+1])/ dx[i]; // Electric field calculation
        }
    }

}

void ES_solver_EC::integrate_time_step(const int thread_id, double del_t, double current_time, const domain& world, std::vector<charged_particle>& particle_list) {
    #pragma omp barrier
    #pragma omp master
    {
        this->particle_timer = MPI_Wtime();
    }
    this->push_particles(thread_id, del_t, particle_list, world);
    this->deposit_charge_density(world, particle_list, thread_id);
    #pragma omp barrier
    #pragma omp master
    {
        double end_time = MPI_Wtime();
        this->particle_timer = end_time - this->particle_timer;
        double start_time = MPI_Wtime();
        this->solve_potential(current_time + del_t, world);
        this->make_EField(world);
        end_time = MPI_Wtime();
        this->potential_timer = end_time - start_time;
    }

}

void ES_solver_EC::push_particles(const int thread_id, double del_t, std::vector<charged_particle>& particle_list, const domain& world){
    // Loop over all particles and push them to the grid
    int num_particles = particle_list.size();
    int left_boundary = world.left_boundary_condition; // Get left boundary condition
    int right_boundary = world.right_boundary_condition; // Get right boundary condition
    int number_cells = world.number_cells; // Number of cells in the domain
    int domain_type = world.domain_type;
    
    if (domain_type == 0) {
        double inv_dx = 1.0 / world.min_dx; // Cell size
        for (int i = 0; i < num_particles; ++i) {
            charged_particle& particle = particle_list[i];
            particle.ES_push_EC_uniform(thread_id, del_t, this->E_field, inv_dx, left_boundary, right_boundary, number_cells); // Push particles to the grid
        }
    } else {
        for (int i = 0; i < num_particles; ++i) {
            charged_particle& particle = particle_list[i];
            particle.ES_push_EC_non_uniform(thread_id, del_t, this->E_field, world.dx_dxi, world.grid_nodes, left_boundary, right_boundary, number_cells); // Push particles to the grid
        }
    }
}



