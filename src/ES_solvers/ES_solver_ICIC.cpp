#include <vector>
#include <omp.h>
#include "ES_solvers/ES_solver_ICIC.hpp"
#include "non_linear_solvers/non_linear_solver.hpp"
#include "globals/mpi_vars.hpp"
#include "globals/constants.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "globals/write_functions.hpp"

ES_solver_ICIC::ES_solver_ICIC(const domain& world) {
    const int number_cells = world.number_cells;
    const int number_nodes = world.number_nodes;
    // potential and charge live on the cell centers, the field on the cell edges
    this->phi.resize(number_cells, 0.0);
    this->phi_past.resize(number_cells, 0.0);
    this->rho.resize(number_cells, 0.0);
    this->E_field.resize(number_nodes, 0.0);
    int number_threads = omp_get_max_threads();
    this->work_space.resize(number_threads);
    for (int i = 0; i < number_threads; i++) {
        this->work_space[i].resize(number_cells, 0.0);
    }
    this->poisson_solver = std::make_unique<poisson_solver_1D_CIC_tridiag>(world);

    // grid metrics used by the field construction and by the mover
    this->dx_cells.resize(number_cells, 0.0);
    for (int i = 0; i < number_cells; ++i) {
        this->dx_cells[i] = (world.domain_type == 0) ? world.min_dx : world.dx_dxi[i];
    }
    this->center_distance.resize(number_cells-1, 0.0);
    for (int i = 0; i < number_cells-1; ++i) {
        this->center_distance[i] = world.cell_centers[i+1] - world.cell_centers[i];
    }
    this->half_dx_left = world.cell_centers[0] - world.grid_nodes[0];
    this->half_dx_right = world.grid_nodes[number_cells] - world.cell_centers[number_cells-1];
    this->wrap_distance = this->half_dx_left + this->half_dx_right;
    // length spanned by each field node's value (see the header note)
    this->node_scale.resize(number_nodes, 1.0);
    for (int i = 1; i < number_cells; ++i) {
        this->node_scale[i] = this->center_distance[i-1];
    }
    if (world.left_boundary_condition == 3) {
        this->node_scale[0] = this->wrap_distance;
        this->node_scale[number_cells] = this->wrap_distance;
    } else {
        // at a wall the field is taken over the half cell, so the difference spans a full cell
        this->node_scale[0] = this->dx_cells[0];
        this->node_scale[number_cells] = this->dx_cells[number_cells-1];
    }
    this->left_boundary_potential = 0.0;
    this->right_boundary_potential = 0.0;
    this->left_boundary_potential_past = 0.0;
    this->right_boundary_potential_past = 0.0;
    this->particle_work_space.resize(number_threads);
    for (int i = 0; i < number_threads; i++) {
        this->particle_work_space[i].resize(number_cells, 0.0);
    }

    std::vector<double> double_params(3);
    std::vector<int> int_params(3);
    read_non_linear_solver_inputs("../inputs/implicit_solver.inp", int_params, double_params);
    MPI_Bcast(double_params.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(int_params.data(), 3, MPI_INT, 0, MPI_COMM_WORLD);
    this->implicit_solver = std::make_unique<AA_solver>(double_params[0], double_params[1], double_params[2], int_params[2], int_params[1], number_cells);
    this->implicit_solver->print_out();

    int flag;
    if (mpi_vars::mpi_rank == 0) {
        std::string line;
        std::ifstream file("../inputs/geometry.inp");
        if (!file) {
            std::cerr << "Error: Unable to open file " << "../inputs/geometry.inp" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for (int i = 0; i < 7; ++i) { std::getline(file, line); }
        std::istringstream iss(line);
        iss >> flag;
        file.close();
    }
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    this->smoothing = (flag == 1);
}

void ES_solver_ICIC::print_out() {
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "ES_solver_ICIC: " << std::endl;
        std::cout << "-------------------------- " << std::endl;
        std::cout << "Number of phi (cell centered) nodes: " << this->phi.size() << std::endl;
        std::cout << "Number of field (cell edge) nodes: " << this->E_field.size() << std::endl;
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

void ES_solver_ICIC::initialize_diagnostic_files(const std::string& filename) {
    if (mpi_vars::mpi_rank == 0) {
        std::ofstream file(filename + "/phi/parameters.dat");
        if (!file) {
            std::cerr << "Error opening file for domain \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        file << "RF_rad_frequency, RF_half_amplitude, left_voltage, right_voltage smoothing \n";
        file << std::scientific << std::setprecision(8);
        file << this->RF_rad_frequency << "\t"
        << this->RF_half_amplitude << "\t"
        << this->left_voltage << "\t"
        << this->right_voltage << "\t"
        << (this->smoothing ? 1 : 0) << "\n";
        file.close();

        file.open(filename + "/field_diagnostics.dat");
        if (!file) {
            std::cerr << "Error opening file for domain \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        file << "field energy (J/m^2),  gauss error \n";
        file.close();
        this->implicit_solver->initialize_diagnostic_files(filename);
    }
}

void ES_solver_ICIC::write_diagnostics(const std::string& dir_name, int diag_number) {
    if (mpi_vars::mpi_rank == 0) {
        write_vector_to_binary_file(this->phi, this->phi.size(), dir_name + "/phi/potential_" + std::to_string(diag_number) + ".dat", 0);
        std::ofstream file(dir_name + "/field_diagnostics.dat", std::ios::app);
        if (!file) {
            std::cerr << "Error opening file for domain \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        file << std::scientific << std::setprecision(8);
        file << this->total_field_energy << "\t"
        << this->gauss_error << "\n";
        file.close();
        this->implicit_solver->write_diagnostics(dir_name);
    }
}

void ES_solver_ICIC::set_boundary_potentials(double current_time, const domain& world) {
    this->left_boundary_potential = (world.left_boundary_condition == 4)
        ? this->RF_half_amplitude * std::sin(this->RF_rad_frequency * current_time) : this->left_voltage;
    this->right_boundary_potential = (world.right_boundary_condition == 4)
        ? this->RF_half_amplitude * std::sin(this->RF_rad_frequency * current_time) : this->right_voltage;
}

// Binomial smoothing of the cell centered charge density.  The image used past a Dirichlet
// wall is anti-symmetric and past a reflective wall symmetric, matching how the deposit
// itself folds the part of the cloud that falls outside the domain.
static void smooth_charge_density(std::vector<double>& rho, const domain& world) {
    const int number_cells = world.number_cells;
    std::vector<double> rho_copy = rho;
    for (int i = 1; i < number_cells - 1; i++) {
        rho[i] = 0.25 * (rho_copy[i-1] + 2.0 * rho_copy[i] + rho_copy[i+1]);
    }
    switch (world.left_boundary_condition) {
        case 1:
        case 4:
            rho[0] = 0.25 * (rho_copy[0] + rho_copy[1]);
            break;
        case 2:
            rho[0] = 0.25 * (3.0 * rho_copy[0] + rho_copy[1]);
            break;
        case 3:
            rho[0] = 0.25 * (rho_copy[number_cells-1] + 2.0 * rho_copy[0] + rho_copy[1]);
            break;
    }
    switch (world.right_boundary_condition) {
        case 1:
        case 4:
            rho[number_cells-1] = 0.25 * (rho_copy[number_cells-1] + rho_copy[number_cells-2]);
            break;
        case 2:
            rho[number_cells-1] = 0.25 * (3.0 * rho_copy[number_cells-1] + rho_copy[number_cells-2]);
            break;
        case 3:
            rho[number_cells-1] = 0.25 * (rho_copy[0] + 2.0 * rho_copy[number_cells-1] + rho_copy[number_cells-2]);
            break;
    }
}

static void smooth_field(std::vector<double>& E_field, const domain& world) {
    const int number_cells = world.number_cells;
    std::vector<double> E_copy = E_field;
    for (int i = 1; i < number_cells; i++) {
        E_field[i] = 0.25 * (E_copy[i-1] + 2.0 * E_copy[i] + E_copy[i+1]);
    }
    switch (world.left_boundary_condition) {
        case 1:
        case 4:
            E_field[0] = 0.5 * (E_copy[0] + E_copy[1]);
            break;
        case 2:
            E_field[0] = 0.0;
            break;
        case 3:
            E_field[0] = 0.25 * (E_copy[number_cells-1] + 2.0 * E_copy[0] + E_copy[1]);
            break;
    }
    switch (world.right_boundary_condition) {
        case 1:
        case 4:
            E_field[number_cells] = 0.5 * (E_copy[number_cells] + E_copy[number_cells-1]);
            break;
        case 2:
            E_field[number_cells] = 0.0;
            break;
        case 3:
            E_field[number_cells] = E_field[0];
            break;
    }
}

void ES_solver_ICIC::deposit_charge_density(const domain& world, std::vector<charged_particle>& particle_list, int thread_id) {
    const int total_thread_count = omp_get_max_threads();
    const int number_cells = world.number_cells;
    const int num_particles = particle_list.size();
    std::vector<double>& part_work_space = this->particle_work_space[thread_id];
    std::vector<double>& local_work_space = this->work_space[thread_id];
    std::fill(local_work_space.begin(), local_work_space.end(), 0.0);
    for (int i = 0; i < num_particles; ++i) {
        const charged_particle& particle = particle_list[i];
        std::fill(part_work_space.begin(), part_work_space.begin() + number_cells, 0.0);
        particle.deposit_particles_quadratic(thread_id, part_work_space,
            world.left_boundary_condition, world.right_boundary_condition, number_cells);
        const double q_times_wp = particle.q_times_wp;
        for (int j = 0; j < number_cells; ++j) {
            local_work_space[j] += part_work_space[j] * q_times_wp;
        }
    }
    #pragma omp barrier
    #pragma omp for
    for (int i = 0; i < number_cells; i++) {
        double sum = 0.0;
        for (int i_thread = 0; i_thread < total_thread_count; i_thread++) {
            sum += this->work_space[i_thread][i];
        }
        this->rho[i] = sum;
    }
    #pragma omp barrier
    #pragma omp master
    {
        MPI_Allreduce(MPI_IN_PLACE, this->rho.data(), number_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (this->smoothing) {
            smooth_charge_density(this->rho, world);
        }
    }
    #pragma omp barrier
}

// Time centred field on the cell edges, in V/m.
void ES_solver_ICIC::make_EField(const domain& world) {
    const int number_cells = world.number_cells;
    for (int i = 1; i < number_cells; ++i) {
        this->E_field[i] = 0.5 * (this->phi[i-1] + this->phi_past[i-1] - this->phi[i] - this->phi_past[i])
            / this->center_distance[i-1];
    }
    switch (world.left_boundary_condition) {
        case 1:
        case 4:
            this->E_field[0] = (0.5 * (this->left_boundary_potential + this->left_boundary_potential_past)
                - 0.5 * (this->phi[0] + this->phi_past[0])) / this->half_dx_left;
            break;
        case 2:
            this->E_field[0] = 0.0;
            break;
        case 3:
            this->E_field[0] = 0.5 * (this->phi[number_cells-1] + this->phi_past[number_cells-1]
                - this->phi[0] - this->phi_past[0]) / this->wrap_distance;
            break;
    }
    switch (world.right_boundary_condition) {
        case 1:
        case 4:
            this->E_field[number_cells] = (0.5 * (this->phi[number_cells-1] + this->phi_past[number_cells-1])
                - 0.5 * (this->right_boundary_potential + this->right_boundary_potential_past)) / this->half_dx_right;
            break;
        case 2:
            this->E_field[number_cells] = 0.0;
            break;
        case 3:
            this->E_field[number_cells] = this->E_field[0];
            break;
    }
    if (this->smoothing) {
        // Smooth the potential difference across each node, not the field.  This is what
        // the Fortran filters, and on a uniform grid the two are the same operator up to a
        // constant.  On a stretched grid only this form keeps the smoothed deposit and the
        // smoothed gather adjoint to one another: the divergence relating them carries no
        // metric in these variables, so the (1 2 1)/4 stencils commute with it exactly,
        // which is what the energy conservation argument rests on.
        for (int i = 0; i <= number_cells; ++i) { this->E_field[i] *= this->node_scale[i]; }
        smooth_field(this->E_field, world);
        for (int i = 0; i <= number_cells; ++i) { this->E_field[i] /= this->node_scale[i]; }
    }
}

void ES_solver_ICIC::solve_potential(double current_time, const domain& world) {
    this->set_boundary_potentials(current_time, world);
    const double inv_epsilon_0 = 1.0 / constants::epsilon_0;
    const int number_cells = world.number_cells;
    for (int i = 0; i < number_cells; ++i) {
        this->phi[i] = -this->rho[i] * inv_epsilon_0;
    }
    // a wall potential enters the first/last row through the half cell gradient
    if (world.left_boundary_condition != 2) {
        this->phi[0] -= this->left_boundary_potential / this->half_dx_left;
    }
    if (world.right_boundary_condition != 2) {
        this->phi[number_cells-1] -= this->right_boundary_potential / this->half_dx_right;
    }
    this->poisson_solver->solve(this->phi, this->phi);
}

void ES_solver_ICIC::solve_field_energy(const domain& world) {
    const int number_cells = world.number_cells;
    double sum = 0.0;
    for (int i = 0; i < number_cells-1; ++i) {
        const double diff = this->phi[i] - this->phi[i+1];
        sum += diff * diff / this->center_distance[i];
    }
    if (world.left_boundary_condition == 3) {
        // single gap across the periodic seam
        const double diff = this->phi[number_cells-1] - this->phi[0];
        sum += diff * diff / this->wrap_distance;
    } else {
        if (world.left_boundary_condition != 2) {
            const double diff = this->left_boundary_potential - this->phi[0];
            sum += diff * diff / this->half_dx_left;
        }
        if (world.right_boundary_condition != 2) {
            const double diff = this->phi[number_cells-1] - this->right_boundary_potential;
            sum += diff * diff / this->half_dx_right;
        }
    }
    this->total_field_energy = 0.5 * constants::epsilon_0 * sum; // J/m^2
}

void ES_solver_ICIC::get_diagnostics(const domain& world, std::vector<charged_particle>& particle_list) {
    const int number_cells = world.number_cells;
    this->solve_field_energy(world);
    // Re-deposit the charge at the final particle positions to check that the converged
    // potential really does satisfy Gauss's law.  This has to use the solver's own deposit
    // rather than the particle density diagnostic: the two differ at a Dirichlet wall,
    // where the solver folds the part of a cloud lying outside the domain back in with the
    // opposite sign as the image charge requires, while the diagnostic folds it in with the
    // same sign so that it still counts every particle.
    const int number_threads = omp_get_max_threads();
    std::vector<double> charge_density(number_cells, 0.0);
    std::vector<double> deposit(number_cells, 0.0);
    for (size_t part_num = 0; part_num < particle_list.size(); part_num++){
        const charged_particle& particle = particle_list[part_num];
        std::fill(deposit.begin(), deposit.end(), 0.0);
        for (int thread_id = 0; thread_id < number_threads; thread_id++) {
            particle.deposit_particles_quadratic(thread_id, deposit,
                world.left_boundary_condition, world.right_boundary_condition, number_cells);
        }
        for (int i = 0; i < number_cells; i++) {
            charge_density[i] += deposit[i] * particle.q_times_wp;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, charge_density.data(), number_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (this->smoothing) {
        smooth_charge_density(charge_density, world);
    }
    std::vector<double> source_term(number_cells, 0.0);
    const double inv_epsilon_0 = 1.0 / constants::epsilon_0;
    for (int i = 0; i < number_cells; ++i) {
        source_term[i] = -charge_density[i] * inv_epsilon_0;
    }
    if (world.left_boundary_condition != 2) {
        source_term[0] -= this->left_boundary_potential / this->half_dx_left;
    }
    if (world.right_boundary_condition != 2) {
        source_term[number_cells-1] -= this->right_boundary_potential / this->half_dx_right;
    }
    this->gauss_error = this->poisson_solver->norm_error(this->phi, source_term);
}

void ES_solver_ICIC::build_mover_coefficients(std::vector<charged_particle>& particle_list, const domain& world) {
    const int number_cells = world.number_cells;
    const size_t num_particles = particle_list.size();
    if (this->mover_cells.size() != num_particles) {
        this->mover_cells.assign(num_particles, std::vector<CIC_cell_coefficients>(number_cells));
    }
    for (size_t i = 0; i < num_particles; ++i) {
        const double q_over_m = particle_list[i].q_over_m;
        std::vector<CIC_cell_coefficients>& cells = this->mover_cells[i];
        // Store the potential difference across each node, not the field.  The mover divides
        // by the cell width to reach logical coordinates, and that width has to cancel against
        // the length the node's value spans, otherwise the work done on the particles no longer
        // matches the change in field energy once the cells stop being equal.
        for (int j = 0; j < number_cells; ++j) {
            cells[j].q_dphi_left = q_over_m * this->E_field[j] * this->node_scale[j];
            cells[j].q_dphi_right = q_over_m * this->E_field[j+1] * this->node_scale[j+1];
            cells[j].inv_dx = 1.0 / this->dx_cells[j];
        }
        // cap the sub-step at a tenth of the local acceleration-gradient time scale so that
        // the Picard iteration inside the cell converges in a few passes.  A cell with a
        // uniform field has no gradient time scale and so is left uncapped.  The cap is a
        // large finite number rather than an infinity because the build enables fast math.
        constexpr double uncapped = 1e30;
        for (int j = 0; j < number_cells; ++j) {
            const double accel_gradient = std::abs(cells[j].q_dphi_left - cells[j].q_dphi_right) * cells[j].inv_dx;
            cells[j].del_tau_min = (accel_gradient > 0.0)
                ? 0.1 * std::sqrt(this->dx_cells[j] / accel_gradient) : uncapped;
        }
    }
}

void ES_solver_ICIC::push_particles(const int thread_id, double del_t, std::vector<charged_particle>& particle_list, const domain& world) {
    const int number_cells = world.number_cells;
    const int left_boundary = world.left_boundary_condition;
    const int right_boundary = world.right_boundary_condition;
    const int num_particles = particle_list.size();
    std::vector<double>& local_work_space = this->work_space[thread_id];
    std::vector<double>& part_work_space = this->particle_work_space[thread_id];
    std::fill(local_work_space.begin(), local_work_space.end(), 0.0);
    for (int i = 0; i < num_particles; ++i) {
        charged_particle& particle = particle_list[i];
        std::fill(part_work_space.begin(), part_work_space.begin() + number_cells, 0.0);
        particle.ES_push_deposit_ICIC(thread_id, del_t, part_work_space,
            this->mover_cells[i], left_boundary, right_boundary, number_cells);
        const double q_times_wp = particle.q_times_wp;
        for (int j = 0; j < number_cells; j++) {
            local_work_space[j] += part_work_space[j] * q_times_wp;
        }
    }
}

void ES_solver_ICIC::write_particle_densities(const std::string file_path, const std::string filename,
    std::vector<charged_particle>& particle_list, const domain& world) const {
    const int number_cells = world.number_cells;
    const int num_particles = particle_list.size();
    for (int i = 0; i < num_particles; ++i) {
        charged_particle& particle = particle_list[i];
        std::vector<double>& density = particle.density;
        // Written unsmoothed: the raw quadratic accumulation on the cell centred nodes divided
        // by the cell size.  Smoothing belongs in post-processing, so that the file always holds
        // what the particles actually deposited and the filter stays a choice made afterwards.
        for (int j = 0; j < number_cells; j++) {
            density[j] = density[j] * particle.weight / this->dx_cells[j];
        }
        if (mpi_vars::mpi_rank == 0) {
            write_vector_to_binary_file(density, number_cells, file_path + "/charged_particles/" + particle.name + "/density/" + filename, 0);
        }
    }
}

void ES_solver_ICIC::integrate_time_step(const int thread_id, double del_t, double current_time,
    const domain& world, std::vector<charged_particle>& particle_list) {

    const int number_cells = world.number_cells;

    auto integral_function = [&](std::vector<double>& res_output) {
        const int total_thread_count = omp_get_max_threads();
        double part_timer_start = 0.0;
        #pragma omp master
        {
            this->set_boundary_potentials(current_time + del_t, world);
            this->make_EField(world);
            this->build_mover_coefficients(particle_list, world);
            part_timer_start = MPI_Wtime();
        }
        #pragma omp barrier
        this->push_particles(thread_id, del_t, particle_list, world);
        #pragma omp barrier
        #pragma omp for
        for (int i = 0; i < number_cells; i++) {
            double sum = 0.0;
            for (int i_thread = 0; i_thread < total_thread_count; i_thread++) {
                sum += this->work_space[i_thread][i];
            }
            this->rho[i] = sum;
        }
        #pragma omp barrier
        #pragma omp master
        {
            MPI_Allreduce(MPI_IN_PLACE, this->rho.data(), number_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if (this->smoothing) {
                smooth_charge_density(this->rho, world);
            }
            double end_time = MPI_Wtime();
            this->particle_timer += (end_time - part_timer_start);
            const double start_time = MPI_Wtime();
            this->solve_potential(current_time + del_t, world);
            end_time = MPI_Wtime();
            this->potential_timer += (end_time - start_time);
        }
        #pragma omp barrier
    };

    #pragma omp master
    {
        this->step_counter++;
        this->phi_past = this->phi;
        this->set_boundary_potentials(current_time, world);
        this->left_boundary_potential_past = this->left_boundary_potential;
        this->right_boundary_potential_past = this->right_boundary_potential;
        this->particle_timer = 0.0;
        this->potential_timer = 0.0;
    }
    #pragma omp barrier
    const size_t failures_before = this->implicit_solver->non_converged_count;
    this->implicit_solver->solve(this->phi, integral_function);
    #pragma omp master
    {
        if (mpi_vars::mpi_rank == 0 && this->implicit_solver->non_converged_count > failures_before) {
            std::cout << "  ^ that non-convergence was at time step " << this->step_counter << std::endl;
        }
    }
    #pragma omp barrier
    #pragma omp master
    {
        this->set_boundary_potentials(current_time + del_t, world);
        this->make_EField(world);
        this->build_mover_coefficients(particle_list, world);
    }
    #pragma omp barrier

    const int num_particles = particle_list.size();
    const int left_boundary = world.left_boundary_condition;
    const int right_boundary = world.right_boundary_condition;
    std::vector<int> number_sub_steps(num_particles, 0);
    for (int i = 0; i < num_particles; ++i) {
        particle_list[i].ES_push_ICIC(thread_id, del_t, number_sub_steps[i],
            this->mover_cells[i], left_boundary, right_boundary, number_cells);
    }
    #pragma omp barrier
}
