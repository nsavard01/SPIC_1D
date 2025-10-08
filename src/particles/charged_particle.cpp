
#include <vector>
#include "particles/charged_particle.hpp"
#include "globals/mpi_vars.hpp"
#include "globals/constants.hpp"
#include "globals/write_functions.hpp"
#include <omp.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <mpi.h>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <dirent.h>
#include "rand_gen/maxwell_generator.hpp"

// static variables for charged_particle class
std::vector<std::vector<double>> charged_particle::xi_sorted, charged_particle::y_sorted, charged_particle::z_sorted, charged_particle::v_x_sorted, charged_particle::v_y_sorted, charged_particle::v_z_sorted;
std::vector<std::vector<size_t>> charged_particle::sorted_number_particles_per_cell, charged_particle::cell_indices; //, charged_particle::particle_cell;

charged_particle::charged_particle(double mass_in, double charge_in, size_t number_in, size_t final_in, std::string name_in, int number_nodes)
    {
    int number_threads = omp_get_max_threads();

    this->name = name_in;
    this->mass = mass_in;
    this->charge = charge_in;
    this->weight = 0.0;
    this->q_over_m = charge_in/mass_in;
    this->q_times_wp = 0.0;
    this->accum_wall_loss[0] = this->accum_wall_loss[1] = 0;
    this->accum_wall_energy_loss[0] = this->accum_wall_energy_loss[0] = 0.0;
    this->total_sum_v[0] = this->total_sum_v[1] = this->total_sum_v[2]  = 0.0;
    this->total_sum_v_square[0] = this->total_sum_v_square[1] = this->total_sum_v_square[2]  = 0.0;
    this->final_idx.resize(number_threads);
    this->number_particles.resize(number_threads);
    this->number_collidable_particles.resize(number_threads);
    this->momentum_loss.resize(number_threads);
    this->accum_wall_momentum_loss.resize(2);
    this->accum_wall_momentum_loss[0].resize(3, 0.0);
    this->accum_wall_momentum_loss[1].resize(3, 0.0);
    this->energy_loss.resize(number_threads);
    this->wall_loss.resize(number_threads);
    this->density.resize(number_nodes, 0.0);
    this->temperature.resize(number_nodes-1, 0.0);
    this->number_particles_per_cell.resize(number_nodes-1, 0);
    if (charged_particle::sorted_number_particles_per_cell.empty()) {
        charged_particle::sorted_number_particles_per_cell.resize(number_threads);
        charged_particle::cell_indices.resize(number_threads);
        // charged_particle::particle_cell.resize(number_threads);
    }
    this->number_unique_injections = 0;
    for (int i = 0; i < number_threads; i++){
        this->number_particles[i].resize(1, number_in);
        this->number_collidable_particles[i].resize(1,number_in);
        this->final_idx[i].resize(1, final_in);
        this->momentum_loss[i].resize(2);
        this->momentum_loss[i][0].resize(3,0.0);
        this->momentum_loss[i][1].resize(3,0.0);
        this->energy_loss[i].resize(2, 0.0);
        this->wall_loss[i].resize(2, 0);
        if (charged_particle::sorted_number_particles_per_cell[i].empty()) {
            charged_particle::sorted_number_particles_per_cell[i].resize(number_nodes-1, 0);
            charged_particle::cell_indices[i].resize(number_nodes, 0);
            // charged_particle::particle_cell[i].resize(final_in, 0);
        }
    }
    
 
}

void charged_particle::print_out() const {
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "Particle name: " << this->name << std::endl;
        std::cout << "Mass (kg): " << this->mass << std::endl;
        std::cout << "Charge (C): " << this->charge << std::endl;
        std::cout << "Weight (m^-2): " << this->weight << std::endl;
        std::cout << "q/m: " << this->q_over_m << std::endl;
        std::cout << "Number of particles per thread: " << this->number_particles[0][0] << std::endl;
        std::cout << "Number of total particles: " << this->total_number_particles << std::endl;
        std::cout << "Final index: " << this->final_idx[0][0] << std::endl;
        std::cout << "Degrees in space: " << this->number_space_coordinates << std::endl;
        std::cout << "Degrees in velocity: " << this->number_velocity_coordinates << std::endl;
        std::cout << "Average density (m^-3): " << this->average_density << std::endl;
        std::cout << "Average temperature (eV): " << this->average_temperature << std::endl;
        std::cout << " " << std::endl;
    }
}

void charged_particle::initialize_number_coordinates(int space, int velocity)
    {
    int number_threads = omp_get_max_threads();
    this->number_space_coordinates = space;
    this->number_velocity_coordinates = velocity;
    
    // add number of spatial and velocity degrees of freedom
    for (int i = 0; i < space; i++){ 
        if (i == 0) {
            this->xi.resize(number_threads);
            for (int j = 0; j < number_threads; j++){
                this->xi[j].resize(this->final_idx[j][0], 0);
            }
            if (charged_particle::xi_sorted.empty()) {
                charged_particle::xi_sorted = this->xi;
            }
        } else if (i == 1) {
            this->y.resize(number_threads);
            for (int j = 0; j < number_threads; j++){
                this->y[j].resize(this->final_idx[j][0], 0);
            }
            if (charged_particle::y_sorted.empty()) {
                charged_particle::y_sorted = this->y;
            }
        } else if (i == 2) {
            this->z.resize(number_threads);
            for (int j = 0; j < number_threads; j++){
                this->z[j].resize(this->final_idx[j][0], 0);
            }
            if (charged_particle::z_sorted.empty()) {
                charged_particle::z_sorted = this->z;
            }
        }
    }
    for (int i = 0; i < velocity; i++){
        if (i == 0) {
            this->v_x.resize(number_threads);
            for (int j = 0; j < number_threads; j++){
                this->v_x[j].resize(this->final_idx[j][0], 0);
            }
            if (charged_particle::v_x_sorted.empty()) {
                charged_particle::v_x_sorted = this->v_x;
            }
        } else if (i == 1) {
            this->v_y.resize(number_threads);
            for (int j = 0; j < number_threads; j++){
                this->v_y[j].resize(this->final_idx[j][0], 0);
            }
            if (charged_particle::v_y_sorted.empty()) {
                charged_particle::v_y_sorted = this->v_y;
            }
        } else if (i == 2) {
            this->v_z.resize(number_threads);
            for (int j = 0; j < number_threads; j++){
                this->v_z[j].resize(this->final_idx[j][0], 0);
            }
            if (charged_particle::v_z_sorted.empty()) {
                charged_particle::v_z_sorted = this->v_z;
            }
        }
    }
    
 
}



void charged_particle::initialize_weight(double n_ave, double L_domain) {
    size_t total_part;
    int i;
    int number_threads = omp_get_max_threads();
    this->average_density = n_ave;
    total_part = 0;
    for (i = 0; i < number_threads; i++){
       total_part += this->number_particles[i][0];
    }
    MPI_Allreduce(&total_part, &this->total_number_particles, 1, mpi_vars::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    this->weight = n_ave * L_domain / static_cast<double>(this->total_number_particles);
    this->q_times_wp = this->charge * this->weight;
}

void charged_particle::initialize_rand_maxwellian(double T_ave, double v_drift) {
    double v_therm = std::sqrt(T_ave * constants::elementary_charge / this->mass);
    double sum_v_sq_x = 0.0;
    double sum_v_sq_y = 0.0;
    double sum_v_sq_z = 0.0;
    double sum_v_x = 0.0;
    double sum_v_y = 0.0;
    double sum_v_z = 0.0;
    if (this->number_velocity_coordinates == 1) {
        #pragma omp parallel reduction(+:sum_v_x, sum_v_sq_x)
        {
            int thread_id = omp_get_thread_num();
            size_t end_part = this->number_particles[thread_id][0];
            for (size_t part_num = 0; part_num < end_part; part_num++){
                double& v_x = this->v_x[thread_id][part_num];
                maxwellian_1D(v_x, v_therm, v_drift);
                sum_v_x += v_x;
                sum_v_sq_x += v_x * v_x;
            }
        }  
        this->total_sum_v[0] = sum_v_x; 
        this->total_sum_v_square[0] = sum_v_sq_x;
    } else if (this->number_velocity_coordinates == 2) {
        #pragma omp parallel reduction(+:sum_v_x, sum_v_y, sum_v_sq_x, sum_v_sq_y)
        {
            int thread_id = omp_get_thread_num();
            size_t end_part = this->number_particles[thread_id][0];
            for (size_t part_num = 0; part_num < end_part; part_num++){
                double& v_x = this->v_x[thread_id][part_num];
                double& v_y = this->v_y[thread_id][part_num];
                maxwellian_2D(v_x, v_y, v_therm, v_drift);
                sum_v_x += v_x;
                sum_v_y += v_y;
                sum_v_sq_x += v_x * v_x; 
                sum_v_sq_y += v_y * v_y;
            }
        }
        this->total_sum_v[0] = sum_v_x; 
        this->total_sum_v[1] = sum_v_y;
        this->total_sum_v_square[0] = sum_v_sq_x;
        this->total_sum_v_square[1] = sum_v_sq_y;

    } else if (this->number_velocity_coordinates == 3) {
        #pragma omp parallel reduction(+:sum_v_sq_x, sum_v_sq_y, sum_v_sq_z, sum_v_x, sum_v_y, sum_v_z)
        {
            int thread_id = omp_get_thread_num();
            size_t end_part = this->number_particles[thread_id][0];
            for (size_t part_num = 0; part_num < end_part; part_num++){
                double& v_x = this->v_x[thread_id][part_num];
                double& v_y = this->v_y[thread_id][part_num];
                double& v_z = this->v_z[thread_id][part_num];
                maxwellian_3D(v_x, v_y, v_z, v_therm, v_drift);
                sum_v_x += v_x;
                sum_v_y += v_y;
                sum_v_z += v_z;
                sum_v_sq_x += v_x * v_x; 
                sum_v_sq_y += v_y * v_y;
                sum_v_sq_z += v_z * v_z;
            }
        }  
        this->total_sum_v[0] = sum_v_x; 
        this->total_sum_v[1] = sum_v_y;
        this->total_sum_v[2] = sum_v_z;
        this->total_sum_v_square[0] = sum_v_sq_x;
        this->total_sum_v_square[1] = sum_v_sq_y;
        this->total_sum_v_square[2] = sum_v_sq_z;
    }
    MPI_Allreduce(MPI_IN_PLACE, this->total_sum_v_square, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->total_sum_v, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double total_sum_v_square_temp = this->total_sum_v_square[0] + this->total_sum_v_square[1] + this->total_sum_v_square[2];
    this->v_sqr_max = 5.0 * (total_sum_v_square_temp) / double(this->total_number_particles); // Assume max about 5x the mean
    this->v_sqr_min = 3.0 * 0.025 / this->mass; // assume minimum at room temperature of 0.025 eV
    this->average_temperature = this->mass * total_sum_v_square_temp / static_cast<double>(this->total_number_particles) / constants::elementary_charge / double(this->number_velocity_coordinates);
}


void charged_particle::initialize_rand_position_uniform(const domain& world) {
    
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        double particles_accum = 0.0;
        double L_domain = world.length_domain;
        size_t particle_total = 0;
        size_t number_part = this->number_particles[thread_id][0];
        int number_cells = world.number_cells;
        size_t start_idx = 0;
        size_t end_idx;
        std::vector<double>& xi_local = this->xi[thread_id];
        for (int i = 0; i < number_cells; i++){
            double dx;
            if (world.domain_type == 0) {
                dx = world.min_dx;
            } else if (world.domain_type == 1) {
                dx = world.dx_dxi[i];
            }
            particles_accum += number_part * dx / L_domain;
            double ideal_particles = particles_accum - particle_total;
            size_t num_particles = std::floor(ideal_particles);
            if (pcg32_random_r() < ideal_particles - num_particles) {
                num_particles += 1;
            }
            end_idx = start_idx + num_particles;
            for (int part_idx = start_idx; part_idx < end_idx; part_idx++){
                xi_local[part_idx] = i + pcg32_random_r();
            }
            particle_total += num_particles;
            start_idx = end_idx;
        }
    }

}

void charged_particle::sort_particle(int thread_id, int number_cells) {

    
    
    bool use_vy = (this->number_velocity_coordinates > 1);
    bool use_vz = (this->number_velocity_coordinates > 2);
    bool use_y = (this->number_space_coordinates > 1);
    bool use_z = (this->number_space_coordinates > 2);
    
    size_t last_idx = this->number_particles[thread_id][0];
    size_t local_indx;
    double xi_temp = 0.0, y_temp = 0.0, z_temp = 0.0, v_x_temp = 0.0, v_y_temp = 0.0, v_z_temp = 0.0;
    std::vector<double>& xi_local = this->xi[thread_id];
    std::vector<double>& v_x_local = this->v_x[thread_id];
    std::vector<double>& v_x_sorted_local = charged_particle::v_x_sorted[thread_id];
    std::vector<size_t>& cell_indices_local = charged_particle::cell_indices[thread_id];
    std::vector<double>& xi_sorted_local = charged_particle::xi_sorted[thread_id];
    std::vector<size_t>& number_part_cell_local = charged_particle::sorted_number_particles_per_cell[thread_id];
    std::vector<double> local_v_sqr(number_cells, 0.0);
    std::fill(number_part_cell_local.begin(), number_part_cell_local.end(), 0);

    
    cell_indices_local[0] = 0;
    for (size_t i = 1; i <= number_cells; i++){
        cell_indices_local[i] = cell_indices_local[i-1] + number_part_cell_local[i-1];
    }
    for (size_t part_num = 0; part_num < last_idx; part_num++){
        xi_temp = xi_local[part_num];
        v_x_temp = v_x_local[part_num];
        local_indx = int(xi_temp);
        size_t pos = cell_indices_local[local_indx]; // current first index of the cell
        xi_sorted_local[pos] = xi_temp;
        v_x_sorted_local[pos] = v_x_temp;
        if (use_y) {
            y_temp = this->y[thread_id][part_num];
            charged_particle::y_sorted[thread_id][pos] = y_temp;
        }
        if (use_z) {
            z_temp = this->z[thread_id][part_num];
            charged_particle::z_sorted[thread_id][pos] = z_temp;
        }
        if (use_vy) {
            v_y_temp = this->v_y[thread_id][part_num];
            charged_particle::v_y_sorted[thread_id][pos] = v_y_temp;
        }
        if (use_vz) {
            v_z_temp = this->v_z[thread_id][part_num];
            charged_particle::v_z_sorted[thread_id][pos] = v_z_temp;
        }
        cell_indices_local[local_indx]++; // increment the first index of the cell
    }
    for (size_t i = 0; i < number_cells; i++){
        cell_indices_local[i] = cell_indices_local[i] + number_part_cell_local[i];
    }
    // swap the sorted arrays with the original ones
    std::swap(xi_local, xi_sorted_local);
    std::swap(v_x_local, v_x_sorted_local);
    if (use_y) {
        std::swap(this->y[thread_id], charged_particle::y_sorted[thread_id]);
    }
    if (use_z) {
        std::swap(this->z[thread_id], charged_particle::z_sorted[thread_id]);
    }
    if (use_vy) {
        std::swap(this->v_y[thread_id], charged_particle::v_y_sorted[thread_id]);
    }
    if (use_vz) {
        std::swap(this->v_z[thread_id], charged_particle::v_z_sorted[thread_id]);
    }

    
}

void charged_particle::get_particle_diagnostics(const int thread_id, const int number_cells, const int density_interp_order) {

    
    // interpolation based on order given
    // for T_e do bin within cell
    #pragma omp master
    {
        this->v_sqr_min = 1e10;
        this->v_sqr_max = 0.0;
    }
    double sum_v_sq[3];
    double sum_v[3];
    sum_v[0] = 0.0; sum_v[1] = 0.0; sum_v[2] = 0.0;
    sum_v_sq[0] = 0.0; sum_v_sq[1] = 0.0; sum_v_sq[2] = 0.0;
    bool use_vy = (this->number_velocity_coordinates > 1);
    bool use_vz = (this->number_velocity_coordinates > 2);
    
    const size_t last_idx = this->number_particles[thread_id][0];
    size_t local_indx;
    double xi_temp = 0.0, v_x_temp = 0.0, v_y_temp = 0.0, v_z_temp = 0.0;
    double v_sqr;
    const std::vector<double>& xi_local = this->xi[thread_id];
    const std::vector<double>& v_x_local = this->v_x[thread_id];
    std::vector<size_t>& number_part_cell_local = charged_particle::sorted_number_particles_per_cell[thread_id];
    std::vector<double> local_density(number_cells+1, 0.0);
    std::vector<double> local_v_sqr(number_cells, 0.0);
    std::fill(number_part_cell_local.begin(), number_part_cell_local.end(), 0);
    double v_sqr_min_local = 1e10;
    double v_sqr_max_local = 0.0;
    for (size_t part_num = 0; part_num < last_idx; part_num++){
        xi_temp = xi_local[part_num];
        v_x_temp = v_x_local[part_num];
        local_indx = int(xi_temp);
        if (density_interp_order == 1) {
            double d = xi_temp - local_indx;
            local_density[local_indx] += (1.0 - d);
            local_density[local_indx+1] += d;
        }
        if (use_vy) {
            v_y_temp = this->v_y[thread_id][part_num];
            sum_v_sq[1] += v_y_temp * v_y_temp;
        }
        if (use_vz) {
            v_z_temp = this->v_z[thread_id][part_num];
            sum_v_sq[2] += v_z_temp * v_z_temp;
        }
        sum_v[0] += v_x_temp;
        sum_v[1] += v_y_temp;
        sum_v[2] += v_z_temp;
        sum_v_sq[0] += v_x_temp * v_x_temp;
        v_sqr = v_x_temp * v_x_temp + v_y_temp * v_y_temp + v_z_temp * v_z_temp;
        v_sqr_min_local = std::min(v_sqr_min_local, v_sqr);
        v_sqr_max_local = std::max(v_sqr_max_local, v_sqr);
        local_v_sqr[local_indx] += v_sqr; // add paticle energy to cell
        number_part_cell_local[local_indx]++;
    }
    


    
    // collect all into net diagnostics
    #pragma omp critical
    {
        this->v_sqr_min = std::min(this->v_sqr_min, v_sqr_min_local);
        this->v_sqr_max = std::max(this->v_sqr_max, v_sqr_max_local);
        for (int i = 0; i < number_cells; i++) {
            this->temperature[i] += local_v_sqr[i];
        }
        for (int i = 0; i < number_cells+1; i++) {
            this->density[i] += local_density[i];
        }
        this->total_sum_v[0] += sum_v[0];
        this->total_sum_v[1] += sum_v[1];
        this->total_sum_v[2] += sum_v[2];
        this->total_sum_v_square[0] += sum_v_sq[0];
        this->total_sum_v_square[1] += sum_v_sq[1];
        this->total_sum_v_square[2] += sum_v_sq[2];
        
    }
    #pragma omp barrier
    #pragma omp for
    for (int i = 0; i < number_cells; i++) {
        for (int i_thread = 0; i_thread < omp_get_max_threads(); i_thread++) {
            this->number_particles_per_cell[i] += charged_particle::sorted_number_particles_per_cell[i_thread][i];
        }
    }
    
    #pragma omp barrier
    
}

void charged_particle::write_diagnostics(const std::string& dir_name, int diag_number) const {
    if (mpi_vars::mpi_rank == 0) {
        
        write_vector_to_binary_file(this->temperature, this->temperature.size(), dir_name + "/charged_particles/" + this->name + "/temperature/cell_temp_" + std::to_string(diag_number) + ".dat", 0);

        std::ofstream file(dir_name + "/charged_particles/" + this->name + "/momentum_diagnostics.dat", std::ios::app);
        if (!file) {
            std::cerr << "Error opening file for momentum particle \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        file << std::scientific << std::setprecision(8);
        file << this->total_sum_v[0] << "\t"
            << this->total_sum_v[1]  << "\t"
            << this->total_sum_v[2] << "\t"
            << this->accum_wall_momentum_loss[0][0] << "\t"
            << this->accum_wall_momentum_loss[0][1] << "\t"
            << this->accum_wall_momentum_loss[0][2] << "\t"
            << this->accum_wall_momentum_loss[1][0] << "\t"
            << this->accum_wall_momentum_loss[1][1] << "\t"
            << this->accum_wall_momentum_loss[1][2]
            <<"\n";

        file.close();

        file.open(dir_name + "/charged_particles/" + this->name + "/energy_diagnostics.dat", std::ios::app);
        if (!file) {
            std::cerr << "Error opening file for energy particle \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        file << std::scientific << std::setprecision(8);
        double sum_v_sq = this->total_sum_v_square[0] + this->total_sum_v_square[1] + this->total_sum_v_square[2]; 
        file << this->total_sum_v_square[0] << "\t"
            << this->total_sum_v_square[1] << "\t"
            << this->total_sum_v_square[2] << "\t"
            << sum_v_sq << "\t"
            << this->accum_wall_energy_loss[0] << "\t"
            << this->accum_wall_energy_loss[1]
            <<"\n";

        file.close();

        file.open(dir_name + "/charged_particles/" + this->name + "/number_diagnostics.dat", std::ios::app);
        if (!file) {
            std::cerr << "Error opening file for number particle \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        file << this->total_number_particles << "\t"
            << this->accum_wall_loss[0] << "\t"
            << this->accum_wall_loss[1]
            <<"\n";

        file.close();
    }
    this->write_phase_space(dir_name);
}

void charged_particle::write_diagnostics_average(const std::string& dir_name) const {
    if (mpi_vars::mpi_rank == 0) {
        
        write_vector_to_binary_file(this->temperature, this->temperature.size(), dir_name + "/charged_particles/" + this->name + "/temperature/cell_temp_average.dat", 0);
        std::ofstream file(dir_name + "/charged_particles/" + this->name + "/momentum_diagnostics_average.dat");
        if (!file) {
            std::cerr << "Error opening file for momentum particle \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        file << "sum_v_x, sum_v_y, sum_v_z, left_sum_v_x, left_sum_v_y, left_sum_v_z, right_sum_v_x, right_sum_v_y, right_sum_v_z \n";
        file << std::scientific << std::setprecision(8);
        file << this->total_sum_v[0] << "\t"
            << this->total_sum_v[1]  << "\t"
            << this->total_sum_v[2] << "\t"
            << this->accum_wall_momentum_loss[0][0] << "\t"
            << this->accum_wall_momentum_loss[0][1] << "\t"
            << this->accum_wall_momentum_loss[0][2] << "\t"
            << this->accum_wall_momentum_loss[1][0] << "\t"
            << this->accum_wall_momentum_loss[1][1] << "\t"
            << this->accum_wall_momentum_loss[1][2]
            <<"\n";

        file.close();

        file.open(dir_name + "/charged_particles/" + this->name + "/energy_diagnostics_average.dat");
        if (!file) {
            std::cerr << "Error opening file for energy particle \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        file << "sum_v_sq_x, sum_v_sq_y, sum_v_sq_z, sum_v_sq, right_sum_v_sq, left_sum_v_sq \n";
        file << std::scientific << std::setprecision(8);
        double sum_v_sq = this->total_sum_v_square[0] + this->total_sum_v_square[1] + this->total_sum_v_square[2]; 
        file << this->total_sum_v_square[0] << "\t"
            << this->total_sum_v_square[1] << "\t"
            << this->total_sum_v_square[2] << "\t"
            << sum_v_sq << "\t"
            << this->accum_wall_energy_loss[0] << "\t"
            << this->accum_wall_energy_loss[1]
            <<"\n";

        file.close();

        file.open(dir_name + "/charged_particles/" + this->name + "/number_diagnostics_average.dat");
        if (!file) {
            std::cerr << "Error opening file for number particle \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        file << "N_p, left_lost, right_lost \n";
        file << this->total_number_particles << "\t"
            << this->accum_wall_loss[0] << "\t"
            << this->accum_wall_loss[1]
            <<"\n";

        file.close();
    }
}

void charged_particle::reset_diagnostics(int thread_id) {
    int density_size = this->density.size();
    int temp_size = this->temperature.size();
    this->wall_loss[thread_id][0] = this->wall_loss[thread_id][1] = 0;
    this->energy_loss[thread_id][0] = this->energy_loss[thread_id][1] = 0;
    for (int j = 0; j < this->number_velocity_coordinates; j++) {
        this->momentum_loss[thread_id][0][j] = 0;
        this->momentum_loss[thread_id][0][j] = 0;
    }
    #pragma omp barrier
    #pragma omp for
    for (int i = 0; i < temp_size; i++){
        this->temperature[i] = 0.0;
        this->number_particles_per_cell[i] = 0;
    }
    #pragma omp for
    for (int i = 0; i < density_size; i++){
        this->density[i] = 0.0;
    }
    #pragma omp master
    {   
        for (int i = 0; i < this->number_velocity_coordinates; i++){
            this->total_sum_v[i] = 0.0;
            this->total_sum_v_square[i] = 0.0;
        }
        this->total_number_particles = 0;
    }
   
}

void charged_particle::write_phase_space(const std::string& dir_name) const {
    // dir_name = file/charged_particles
    std::string file_path = dir_name + "/charged_particles/" + this->name + "/phase_space/";
    for (int rank_num = 0; rank_num < mpi_vars::mpi_size; rank_num++){
        if (mpi_vars::mpi_rank == rank_num) {
            for (int i_thread = 0; i_thread<omp_get_max_threads();i_thread++){
                bool append = !(mpi_vars::mpi_rank == 0 && i_thread == 0); // initialize for rank =0, thread = 0
                size_t part_num = this->number_particles[i_thread][0];
                write_vector_to_binary_file(this->xi[i_thread], part_num, file_path + "xi.dat", 0, 0x00, append);
                write_vector_to_binary_file(this->v_x[i_thread], part_num, file_path + "v_x.dat", 0, 0x00, append);
                if (this->number_space_coordinates > 1) {
                    write_vector_to_binary_file(this->y[i_thread], part_num, file_path + "y.dat", 0, 0x00, append);
                }
                if (this->number_space_coordinates > 2) {
                    write_vector_to_binary_file(this->z[i_thread], part_num, file_path + "z.dat", 0, 0x00, append);
                }
                if (this->number_velocity_coordinates > 1) {
                    write_vector_to_binary_file(this->v_y[i_thread], part_num, file_path + "v_y.dat", 0, 0x00, append);
                }
                if (this->number_velocity_coordinates > 2) {
                    write_vector_to_binary_file(this->v_z[i_thread], part_num, file_path + "v_z.dat", 0, 0x00, append);
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void charged_particle::gather_mpi(){
    // Should be done after sorting and diagnostics
    
    this->total_number_particles = 0;
    this->accum_wall_loss[0] = this->accum_wall_loss[1] = 0;
    this->accum_wall_energy_loss[0] = this->accum_wall_energy_loss[1] = 0;
    for (int i = 0; i < this->number_velocity_coordinates; i++) {
        this->accum_wall_momentum_loss[0][i] = this->accum_wall_momentum_loss[1][i] = 0;
    }
    for (int i = 0; i < omp_get_max_threads(); i++){
        this->total_number_particles += this->number_particles[i][0];
        this->accum_wall_loss[0] += this->wall_loss[i][0];
        this->accum_wall_loss[1] += this->wall_loss[i][1];
        this->accum_wall_energy_loss[0] += this->energy_loss[i][0];
        this->accum_wall_energy_loss[1] += this->energy_loss[i][1];
        for (int j = 0; j < this->number_velocity_coordinates; j++) {
            this->accum_wall_momentum_loss[0][j] += this->momentum_loss[i][0][j];
            this->accum_wall_momentum_loss[1][j] += this->momentum_loss[i][1][j];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &this->total_number_particles, 1, mpi_vars::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->total_sum_v, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->total_sum_v_square, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &this->v_sqr_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &this->v_sqr_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->accum_wall_energy_loss, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->accum_wall_momentum_loss[0].data(), 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->accum_wall_momentum_loss[1].data(), 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->accum_wall_loss, 2, mpi_vars::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->temperature.data(), this->temperature.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->density.data(), this->density.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(MPI_IN_PLACE, this->number_particles_per_cell.data(), this->number_particles_per_cell.size(), mpi_vars::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    this->average_temperature = 0.0;
    size_t total_number_particles_local = 0;
    for (size_t i = 0; i < this->temperature.size(); i++) {
        this->average_temperature += this->temperature[i];
        total_number_particles_local += this->number_particles_per_cell[i];
        this->temperature[i] = this->temperature[i] * this->mass /static_cast<double>(this->number_particles_per_cell[i])/ double(this->number_velocity_coordinates) / constants::elementary_charge; // convert to temperature     
    }
    this->average_temperature = this->average_temperature * this->mass / static_cast<double>(total_number_particles_local) / constants::elementary_charge / double(this->number_velocity_coordinates);
    this->average_density = static_cast<double>(total_number_particles_local);
    this->average_density *= this->weight; // convert to density
    
}


void charged_particle::deposit_particles_linear(const int thread_id, std::vector<double>& work_space) const {
    size_t last_idx = this->number_particles[thread_id][0];
    double d, xi_p;
    int xi_left, xi_right;
    const std::vector<double>& xi_local = this->xi[thread_id];
    for (size_t part_indx = 0; part_indx < last_idx; part_indx++){
        xi_p = xi_local[part_indx];
        xi_left = int(xi_p);
        xi_right = xi_left+1;
        d = xi_p - xi_left;
        work_space[xi_left] += (1.0 - d);
        work_space[xi_right] += d;
    } 
}


void charged_particle::ES_push_MC(const int thread_id, double del_t, const std::vector<double>& E_field, const double inv_dx, const int left_boundary, const int right_boundary, const int number_cells) {
   
    size_t last_idx = this->number_particles[thread_id][0];
    std::vector<double>& xi_local = this->xi[thread_id];
    std::vector<double>& v_x_local = this->v_x[thread_id];
    size_t space_delete = 0;
    double v_x, xi, v_y = 0.0, v_z = 0.0;
    double E_field_local, d_p;
    bool use_vy = (this->number_velocity_coordinates > 1);
    bool use_vz = (this->number_velocity_coordinates > 2);
    bool use_y = (this->number_space_coordinates > 1);
    bool use_z = (this->number_space_coordinates > 2);
    int xi_cell;
    bool del_part;
    
    for (size_t part_indx= 0; part_indx < last_idx; part_indx++){
        xi = xi_local[part_indx];
        xi_cell = int(xi);
        d_p = xi - xi_cell;
        E_field_local = E_field[xi_cell] * (1.0 - d_p) + E_field[xi_cell+1] * d_p;
        v_x = v_x_local[part_indx];
        v_x += this->q_over_m * E_field_local * del_t;
        xi += v_x * del_t * inv_dx;
        del_part = false;
        if (use_vy) {
            v_y = this->v_y[thread_id][part_indx];
        }
        if (use_vz){
            v_z = this->v_z[thread_id][part_indx];
        }
        if (xi <= 0) {
            switch (left_boundary){
                case 1:
                case 4:
                    this->energy_loss[thread_id][0] += v_x*v_x + v_y*v_y + v_z*v_z;
                    this->wall_loss[thread_id][0]++;
                    this->momentum_loss[thread_id][0][0] += v_x;
                    this->momentum_loss[thread_id][0][1] += v_y;
                    this->momentum_loss[thread_id][0][2] += v_z;
                    del_part = true;
                    break;
                case 2:
                    xi = - xi;
                    v_x = - v_x;
                    break;
                case 3:
                    xi = number_cells + xi;
                    break;
            }
        } else if (xi >= number_cells) {
            switch (right_boundary){
                case 1:
                case 4:
                    this->energy_loss[thread_id][1] += v_x*v_x + v_y*v_y + v_z*v_z;
                    this->wall_loss[thread_id][1]++;
                    this->momentum_loss[thread_id][1][0] += v_x;
                    this->momentum_loss[thread_id][1][1] += v_y;
                    this->momentum_loss[thread_id][1][2] += v_z;
                    del_part = true;
                    break;
                case 2:
                    xi = 2.0 * number_cells - xi;
                    v_x = - v_x;
                    break;
                case 3:
                    xi = xi - number_cells;
                    break;
            }

        }
        if (!del_part) {
            size_t new_idx = part_indx-space_delete;
            xi_local[new_idx] = xi;
            v_x_local[new_idx] = v_x;
            if (use_vy) {
                this->v_y[thread_id][new_idx] = v_y;
            }
            if (use_vz){
                this->v_z[thread_id][new_idx] = v_z;
            }
            if (use_y) {
                this->y[thread_id][new_idx] = this->y[thread_id][part_indx];
            }
            if (use_z) {
                this->z[thread_id][new_idx] = this->z[thread_id][part_indx];
            }
        } else {
            space_delete++;
        }
    }

    // now for injected particles with different del_t
    for (int inj_indx = 0; inj_indx < this->number_unique_injections; inj_indx++) {
        size_t start_indx = last_idx;
        const size_t number_inject_local = this->number_particles_injected[thread_id][inj_indx];
        const std::vector<double>& del_t_array = this->time_step_injected[thread_id][inj_indx];
        last_idx = start_indx + number_inject_local;
        double del_t_local;
        for (size_t part_indx = start_indx; part_indx < last_idx; part_indx++){
            xi = xi_local[part_indx];
            del_t_local = del_t_array[part_indx-start_indx];
            xi_cell = int(xi);
            d_p = xi - xi_cell;
            E_field_local = E_field[xi_cell] * (1.0 - d_p) + E_field[xi_cell+1] * d_p;
            v_x = v_x_local[part_indx];
            v_x += this->q_over_m * E_field_local * del_t_local;
            xi += v_x * del_t_local * inv_dx;
            del_part = false;
            if (use_vy) {
                v_y = this->v_y[thread_id][part_indx];
            }
            if (use_vz){
                v_z = this->v_z[thread_id][part_indx];
            }
            if (xi <= 0) {
                switch (left_boundary){
                    case 1:
                    case 4:
                        this->energy_loss[thread_id][0] += v_x*v_x + v_y*v_y + v_z*v_z;
                        this->wall_loss[thread_id][0]++;
                        this->momentum_loss[thread_id][0][0] += v_x;
                        this->momentum_loss[thread_id][0][1] += v_y;
                        this->momentum_loss[thread_id][0][2] += v_z;
                        del_part = true;
                        break;
                    case 2:
                        xi = - xi;
                        v_x = - v_x;
                        break;
                    case 3:
                        xi = number_cells + xi;
                        break;
                }
            } else if (xi >= number_cells) {
                switch (right_boundary){
                    case 1:
                    case 4:
                        this->energy_loss[thread_id][1] += v_x*v_x + v_y*v_y + v_z*v_z;
                        this->wall_loss[thread_id][1]++;
                        this->momentum_loss[thread_id][1][0] += v_x;
                        this->momentum_loss[thread_id][1][1] += v_y;
                        this->momentum_loss[thread_id][1][2] += v_z;
                        del_part = true;
                        break;
                    case 2:
                        xi = 2.0 * number_cells - xi;
                        v_x = - v_x;
                        break;
                    case 3:
                        xi = xi - number_cells;
                        break;
                }

            }
            if (!del_part) {
                size_t new_idx = part_indx-space_delete;
                xi_local[new_idx] = xi;
                v_x_local[new_idx] = v_x;
                if (use_vy) {
                    this->v_y[thread_id][new_idx] = v_y;
                }
                if (use_vz){
                    this->v_z[thread_id][new_idx] = v_z;
                }
                if (use_y) {
                    this->y[thread_id][new_idx] = this->y[thread_id][part_indx];
                }
                if (use_z) {
                    this->z[thread_id][new_idx] = this->z[thread_id][part_indx];
                }
            } else {
                space_delete++;
            }
        }
    }
    this->number_particles[thread_id][0] = (last_idx - space_delete);
    this->number_collidable_particles[thread_id][0] = (last_idx - space_delete);
}

void charged_particle::ES_push_EC_uniform(const int thread_id, double del_t, const std::vector<double>& E_field, const double inv_dx, const int left_boundary, const int right_boundary, const int number_cells) {
   
    size_t last_idx = this->number_particles[thread_id][0];
    std::vector<double>& xi_local = this->xi[thread_id];
    std::vector<double>& v_x_local = this->v_x[thread_id];
    size_t space_delete = 0;
    double v_x, xi, v_y = 0.0, v_z = 0.0;
    double E_field_local;
    bool use_vy = (this->number_velocity_coordinates > 1);
    bool use_vz = (this->number_velocity_coordinates > 2);
    bool use_y = (this->number_space_coordinates > 1);
    bool use_z = (this->number_space_coordinates > 2);
    int xi_cell;
    bool del_part;
    for (size_t part_indx= 0; part_indx < last_idx; part_indx++){
        xi = xi_local[part_indx];
        xi_cell = int(xi);
        E_field_local = E_field[xi_cell];
        v_x = v_x_local[part_indx];
        v_x += this->q_over_m * E_field_local * del_t;
        xi += v_x * del_t * inv_dx;
        del_part = false;
        if (use_vy) {
            v_y = this->v_y[thread_id][part_indx];
        }
        if (use_vz){
            v_z = this->v_z[thread_id][part_indx];
        }
        if (xi <= 0) {
            switch (left_boundary){
                case 1:
                case 4:
                    this->energy_loss[thread_id][0] += v_x*v_x + v_y*v_y + v_z*v_z;
                    this->wall_loss[thread_id][0]++;
                    this->momentum_loss[thread_id][0][0] += v_x;
                    this->momentum_loss[thread_id][0][1] += v_y;
                    this->momentum_loss[thread_id][0][2] += v_z;
                    del_part = true;
                    break;
                case 2:
                    xi = - xi;
                    v_x = - v_x;
                    break;
                case 3:
                    xi = number_cells + xi;
                    break;
            }
        } else if (xi >= number_cells) {
            switch (right_boundary){
                case 1:
                case 4:
                    this->energy_loss[thread_id][1] += v_x*v_x + v_y*v_y + v_z*v_z;
                    this->wall_loss[thread_id][1]++;
                    this->momentum_loss[thread_id][1][0] += v_x;
                    this->momentum_loss[thread_id][1][1] += v_y;
                    this->momentum_loss[thread_id][1][2] += v_z;
                    del_part = true;
                    break;
                case 2:
                    xi = 2.0 * number_cells - xi;
                    v_x = - v_x;
                    break;
                case 3:
                    xi = xi - number_cells;
                    break;
            }

        }
        if (!del_part) {
            size_t new_idx = part_indx-space_delete;
            xi_local[new_idx] = xi;
            v_x_local[new_idx] = v_x;
            if (use_vy) {
                this->v_y[thread_id][new_idx] = v_y;
            }
            if (use_vz){
                this->v_z[thread_id][new_idx] = v_z;
            }
            if (use_y) {
                this->y[thread_id][new_idx] = this->y[thread_id][part_indx];
            }
            if (use_z) {
                this->z[thread_id][new_idx] = this->z[thread_id][part_indx];
            }
        } else {
            space_delete++;
        }
    }
    for (int inj_indx = 0; inj_indx < this->number_unique_injections; inj_indx++) {
        size_t start_indx = last_idx;
        const size_t number_inject_local = this->number_particles_injected[thread_id][inj_indx];
        const std::vector<double>& del_t_array = this->time_step_injected[thread_id][inj_indx];
        last_idx = start_indx + number_inject_local;
        double del_t_local;
        for (size_t part_indx = start_indx; part_indx < last_idx; part_indx++){
            xi = xi_local[part_indx];
            del_t_local = del_t_array[part_indx-start_indx];
            xi_cell = int(xi);
            E_field_local = E_field[xi_cell];
            v_x = v_x_local[part_indx];
            v_x += this->q_over_m * E_field_local * del_t_local;
            xi += v_x * del_t_local * inv_dx;
            del_part = false;
            if (use_vy) {
                v_y = this->v_y[thread_id][part_indx];
            }
            if (use_vz){
                v_z = this->v_z[thread_id][part_indx];
            }
            if (xi <= 0) {
                switch (left_boundary){
                    case 1:
                    case 4:
                        this->energy_loss[thread_id][0] += v_x*v_x + v_y*v_y + v_z*v_z;
                        this->wall_loss[thread_id][0]++;
                        this->momentum_loss[thread_id][0][0] += v_x;
                        this->momentum_loss[thread_id][0][1] += v_y;
                        this->momentum_loss[thread_id][0][2] += v_z;
                        del_part = true;
                        break;
                    case 2:
                        xi = - xi;
                        v_x = - v_x;
                        break;
                    case 3:
                        xi = number_cells + xi;
                        break;
                }
            } else if (xi >= number_cells) {
                switch (right_boundary){
                    case 1:
                    case 4:
                        this->energy_loss[thread_id][1] += v_x*v_x + v_y*v_y + v_z*v_z;
                        this->wall_loss[thread_id][1]++;
                        this->momentum_loss[thread_id][1][0] += v_x;
                        this->momentum_loss[thread_id][1][1] += v_y;
                        this->momentum_loss[thread_id][1][2] += v_z;
                        del_part = true;
                        break;
                    case 2:
                        xi = 2.0 * number_cells - xi;
                        v_x = - v_x;
                        break;
                    case 3:
                        xi = xi - number_cells;
                        break;
                }

            }
            if (!del_part) {
                size_t new_idx = part_indx-space_delete;
                xi_local[new_idx] = xi;
                v_x_local[new_idx] = v_x;
                if (use_vy) {
                    this->v_y[thread_id][new_idx] = v_y;
                }
                if (use_vz){
                    this->v_z[thread_id][new_idx] = v_z;
                }
                if (use_y) {
                    this->y[thread_id][new_idx] = this->y[thread_id][part_indx];
                }
                if (use_z) {
                    this->z[thread_id][new_idx] = this->z[thread_id][part_indx];
                }
            } else {
                space_delete++;
            }
        }
    }
    this->number_particles[thread_id][0] = (last_idx - space_delete);
    this->number_collidable_particles[thread_id][0] = (last_idx - space_delete);
}

void charged_particle::ES_push_EC_non_uniform(const int thread_id, double del_t, const std::vector<double>& E_field, 
    const std::vector<double>& dx_dxi, const std::vector<double>& grid, const int left_boundary, const int right_boundary, const int number_cells) {
   
    size_t last_idx = this->number_particles[thread_id][0];
    std::vector<double>& xi_local = this->xi[thread_id];
    std::vector<double>& v_x_local = this->v_x[thread_id];
    size_t space_delete = 0;
    double v_x, xi, x_i, x_f, x_left, v_y = 0.0, v_z = 0.0;
    double E_field_local;
    bool use_vy = (this->number_velocity_coordinates > 1);
    bool use_vz = (this->number_velocity_coordinates > 2);
    bool use_y = (this->number_space_coordinates > 1);
    bool use_z = (this->number_space_coordinates > 2);
    int xi_cell;
    bool del_part;
    double dx;
    const double x_left_boundary = grid[0];
    const double x_right_boundary = grid[number_cells];
    for (size_t part_indx= 0; part_indx < last_idx; part_indx++){
        xi = xi_local[part_indx];
        xi_cell = int(xi);
        E_field_local = E_field[xi_cell];
        v_x = v_x_local[part_indx];
        v_x += this->q_over_m * E_field_local * del_t;
        dx = dx_dxi[xi_cell];
        x_left = grid[xi_cell];
        x_i = x_left + dx * (xi - xi_cell);
        x_f = x_i +  v_x * del_t;
        del_part = false;
        if (use_vy) {
            v_y = this->v_y[thread_id][part_indx];
        }
        if (use_vz){
            v_z = this->v_z[thread_id][part_indx];
        }
        if (x_f <= x_left_boundary) {
            switch (left_boundary){
                case 1:
                case 4:
                    this->energy_loss[thread_id][0] += v_x*v_x + v_y*v_y + v_z*v_z;
                    this->wall_loss[thread_id][0]++;
                    this->momentum_loss[thread_id][0][0] += v_x;
                    this->momentum_loss[thread_id][0][1] += v_y;
                    this->momentum_loss[thread_id][0][2] += v_z;
                    del_part = true;
                    break;
                case 2:
                    x_f = -x_f;
                    v_x = - v_x;
                    break;
                case 3:
                    x_f = x_right_boundary + x_f;
                    break;
            }
        } else if (x_f >= x_right_boundary) {
            switch (right_boundary){
                case 1:
                case 4:
                    this->energy_loss[thread_id][1] += v_x*v_x + v_y*v_y + v_z*v_z;
                    this->wall_loss[thread_id][1]++;
                    this->momentum_loss[thread_id][1][0] += v_x;
                    this->momentum_loss[thread_id][1][1] += v_y;
                    this->momentum_loss[thread_id][1][2] += v_z;
                    del_part = true;
                    break;
                case 2:
                    x_f = 2.0 * x_right_boundary - x_f;
                    v_x = - v_x;
                    break;
                case 3:
                    x_f = x_f - x_right_boundary;
                    break;
            }

        }
        if (!del_part) {
            size_t new_idx = part_indx-space_delete;
            double diff = x_f - x_i;
            int direction = (diff > 0) - (diff < 0);
            while (x_f < grid[xi_cell] || x_f >= grid[xi_cell+1]){
                xi_cell += direction;
            }
            xi_local[new_idx] = xi_cell + (x_f - grid[xi_cell])/dx_dxi[xi_cell];
            v_x_local[new_idx] = v_x;
            if (use_vy) {
                this->v_y[thread_id][new_idx] = v_y;
            }
            if (use_vz){
                this->v_z[thread_id][new_idx] = v_z;
            }
            if (use_y) {
                this->y[thread_id][new_idx] = this->y[thread_id][part_indx];
            }
            if (use_z) {
                this->z[thread_id][new_idx] = this->z[thread_id][part_indx];
            }
        } else {
            space_delete++;
        }
    }
    for (int inj_indx = 0; inj_indx < this->number_unique_injections; inj_indx++) {
        size_t start_indx = last_idx;
        const size_t number_inject_local = this->number_particles_injected[thread_id][inj_indx];
        const std::vector<double>& del_t_array = this->time_step_injected[thread_id][inj_indx];
        last_idx = start_indx + number_inject_local;
        double del_t_local;
        for (size_t part_indx = start_indx; part_indx < last_idx; part_indx++){
            xi = xi_local[part_indx];
            del_t_local = del_t_array[part_indx-start_indx];
            xi_cell = int(xi);
            E_field_local = E_field[xi_cell];
            v_x = v_x_local[part_indx];
            v_x += this->q_over_m * E_field_local * del_t_local;
            dx = dx_dxi[xi_cell];
            x_left = grid[xi_cell];
            x_i = x_left + dx * (xi - xi_cell);
            x_f = x_i +  v_x * del_t_local;
            del_part = false;
            if (use_vy) {
                v_y = this->v_y[thread_id][part_indx];
            }
            if (use_vz){
                v_z = this->v_z[thread_id][part_indx];
            }
            if (x_f <= x_left_boundary) {
                switch (left_boundary){
                    case 1:
                    case 4:
                        this->energy_loss[thread_id][0] += v_x*v_x + v_y*v_y + v_z*v_z;
                        this->wall_loss[thread_id][0]++;
                        this->momentum_loss[thread_id][0][0] += v_x;
                        this->momentum_loss[thread_id][0][1] += v_y;
                        this->momentum_loss[thread_id][0][2] += v_z;
                        del_part = true;
                        break;
                    case 2:
                        x_f = -x_f;
                        v_x = - v_x;
                        break;
                    case 3:
                        x_f = x_right_boundary + x_f;
                        break;
                }
            } else if (x_f >= x_right_boundary) {
                switch (right_boundary){
                    case 1:
                    case 4:
                        this->energy_loss[thread_id][1] += v_x*v_x + v_y*v_y + v_z*v_z;
                        this->wall_loss[thread_id][1]++;
                        this->momentum_loss[thread_id][1][0] += v_x;
                        this->momentum_loss[thread_id][1][1] += v_y;
                        this->momentum_loss[thread_id][1][2] += v_z;
                        del_part = true;
                        break;
                    case 2:
                        x_f = 2.0 * x_right_boundary - x_f;
                        v_x = - v_x;
                        break;
                    case 3:
                        x_f = x_f - x_right_boundary;
                        break;
                }

            }
            if (!del_part) {
                size_t new_idx = part_indx-space_delete;
                double diff = x_f - x_i;
                int direction = (diff > 0) - (diff < 0);
                while (x_f < grid[xi_cell] || x_f >= grid[xi_cell+1]){
                    xi_cell += direction;
                }
                xi_local[new_idx] = xi_cell + (x_f - grid[xi_cell])/dx_dxi[xi_cell];
                v_x_local[new_idx] = v_x;
                if (use_vy) {
                    this->v_y[thread_id][new_idx] = v_y;
                }
                if (use_vz){
                    this->v_z[thread_id][new_idx] = v_z;
                }
                if (use_y) {
                    this->y[thread_id][new_idx] = this->y[thread_id][part_indx];
                }
                if (use_z) {
                    this->z[thread_id][new_idx] = this->z[thread_id][part_indx];
                }
            } else {
                space_delete++;
            }
        }
    }
    this->number_particles[thread_id][0] = (last_idx - space_delete);
    this->number_collidable_particles[thread_id][0] = (last_idx - space_delete);
}



// void Particle::load_density(bool reset_bool){
//     if (reset_bool) {
//         #pragma omp parallel for
//         for (size_t i = 0; i < global_inputs::number_nodes;i++) {
//             this->density[i] = 0.0;
//         }
//     }
//     #pragma omp parallel
//     {
//         this->interpolate_particles();
//         #pragma omp barrier
//         int thread_id = omp_get_thread_num();
//         #pragma omp critical
//         {   
//             for (int i = 0; i<global_inputs::number_nodes; i++) {
//                 this->density[i] += this->work_space[thread_id][i];
//             }
//         }
//     }

// }

// void Particle::write_phase_space(const std::string& dir_name, int diag_num) const{
//     for (int rank_num = 0; rank_num < Constants::mpi_size; rank_num++){
//         if (Constants::mpi_rank == rank_num) {
//             for (int i_thread = 0; i_thread<global_inputs::number_omp_threads;i_thread++){
//                 size_t size_vec = this->number_particles[i_thread]*4;
//                 std::vector<double> temp_phase_space(this->phase_space[i_thread].begin(), this->phase_space[i_thread].begin() + size_vec);
//                 bool append = (Constants::mpi_rank != 0);
                
//                 std::string filename = dir_name + "/PhaseSpace/phaseSpace_" + this->name + "_thread" + std::to_string(i_thread+1) + ".dat";
//                 write_vector_to_binary_file(temp_phase_space, filename, 0, 0x00, append);
//             }
//         }
//         MPI_Barrier(MPI_COMM_WORLD);
//     }
// }

// void Particle::write_density(const std::string& dir_name, const Domain& world, size_t current_diag, bool average_bool){
//     if (world.boundary_conditions[0] == 3) {
//         this->density[0] = this->density[0] + this->density[global_inputs::number_cells];
//         this->density[global_inputs::number_cells] = this->density[0];
//     } else {
//         this->density[0] = 2.0*this->density[0];
//         this->density[global_inputs::number_cells] = 2.0 * this->density[global_inputs::number_cells];
//     }
//     for (int i = 0;i<global_inputs::number_nodes; i++) {
//         this->density[i] = this->density[i] * this->weight/ world.del_X;
//     }
//     MPI_Allreduce(MPI_IN_PLACE, this->density.data(), global_inputs::number_nodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//     std::string filename;
//     if (average_bool) {
//         filename = dir_name + "/Density/density_" + this->name + "_Average.dat";
//         for (int i = 0; i<global_inputs::number_nodes;i++){
//             this->density[i] /= current_diag; //current diag is averaging number if averaging
//         }
//     } else {
//         filename = dir_name + "/Density/density_" + this->name + "_" + std::to_string(current_diag) + ".dat";
//     }
//     if (Constants::mpi_rank == 0) {write_vector_to_binary_file(this->density, filename, 4);}
// }


// double Particle::get_KE_ave() const{
//     // in eV
//     return this->total_sum_v_square * this->mass * 0.5 / this->total_number_particles / Constants::elementary_charge;
// }

// double Particle::get_KE_total() const{
//     return this->total_sum_v_square * this->mass * 0.5 * this->weight;
// }

// double Particle::get_momentum_total() const{
//     return this->total_sum_v_x * this->mass;
// }

// void Particle::write_cell_temperature(const std::string& dir_name, int diag_num) const {
//     std::vector<double> temp(global_inputs::number_cells, 0.0);
//     std::vector<size_t> counter(global_inputs::number_cells, 0);
//     #pragma omp parallel
//     {
//         int thread_id = omp_get_thread_num();
//         std::vector<double> temp_thread(global_inputs::number_cells, 0.0);
//         std::vector<size_t> counter_thread(global_inputs::number_cells, 0);
//         size_t last_idx = this->number_particles[thread_id]*4;
//         double v_x, v_y, v_z;
//         size_t xi_cell;
//         for (size_t part_idx = 0; part_idx < last_idx; part_idx += 4){
//             xi_cell = int(this->phase_space[thread_id][part_idx]);
//             v_x = this->phase_space[thread_id][part_idx + 1];
//             v_y = this->phase_space[thread_id][part_idx + 2];
//             v_z = this->phase_space[thread_id][part_idx + 3];
//             temp_thread[xi_cell] += v_x*v_x + v_y*v_y + v_z*v_z;
//             counter_thread[xi_cell]++;
//         }
//         #pragma omp critical
//         {
//             for (int i = 0; i<global_inputs::number_cells;i++){
//                 temp[i] += temp_thread[i];
//                 counter[i] += counter_thread[i];
//             }
//         }
//     }
//     MPI_Allreduce(MPI_IN_PLACE, temp.data(), global_inputs::number_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//     MPI_Allreduce(MPI_IN_PLACE, counter.data(), global_inputs::number_cells, Constants::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
//     for (int i = 0; i<global_inputs::number_cells;i++){
//         if (counter[i] > 0) {
//             temp[i] = temp[i] * this->mass / counter[i] / 3.0/Constants::elementary_charge;
//         } else {
//             temp[i] = 0.0;
//         }
//     }
//     std::string filename = dir_name + "/Temperature/Temp_" + this->name + "_" + std::to_string(diag_num) + ".dat";
//     if (Constants::mpi_rank == 0) {write_vector_to_binary_file<double>(temp, filename, 4);}

// }



void charged_particle::initialize_diagnostic_files(const std::string& dir_name) const {
    if (mpi_vars::mpi_rank == 0) {

        std::ofstream file(dir_name + "/charged_particles/" + this->name + "/particle_properties.dat");

        // Write header (optional)
        file << "Particle Symbol, Particle Mass (kg), Particle Charge (C), Particle Weight (N/m^2), number space coord, number velocity coord \n";
        file << std::scientific << std::setprecision(8);
        file << this->name << "\t"
            << this->mass << "\t"
            << this->charge << "\t"
            << this->weight << "\t"
            << this->number_space_coordinates << "\t"
            << this->number_velocity_coordinates
            <<"\n";


        file.close();

        file.open(dir_name + "/charged_particles/" + this->name + "/momentum_diagnostics.dat");
        if (!file) {
            std::cerr << "Error opening file for momentum particle \n";
            return;
        }

        file << "sum_v_x, sum_v_y, sum_v_z, left_sum_v_x, left_sum_v_y, left_sum_v_z, right_sum_v_x, right_sum_v_y, right_sum_v_z \n";

        file.close();

        file.open(dir_name + "/charged_particles/" + this->name + "/energy_diagnostics.dat");
        if (!file) {
            std::cerr << "Error opening file for energy particle \n";
            return;
        }

        file << "sum_v_sq_x, sum_v_sq_y, sum_v_sq_z, sum_v_sq, right_sum_v_sq, left_sum_v_sq \n";

        file.close();

        file.open(dir_name + "/charged_particles/" + this->name + "/number_diagnostics.dat");
        if (!file) {
            std::cerr << "Error opening file for energy particle \n";
            return;
        }

        file << "N_p, left_lost, right_lost \n";

        file.close();
    }

}

// void Particle::diag_write(const std::string& dir_name, const double& time_diff, const double& current_time, bool average_bool) const {
//     double KE_ave = this->get_KE_ave() * 2.0/3.0;
//     if (Constants::mpi_rank == 0) {
//         if (!average_bool) {
//             std::ofstream file(dir_name + "/ParticleDiagnostic_" + this->name + ".dat", std::ios::app);
//             if (!file) {
//                 std::cerr << "Error opening file\n";
//                 return;
//             }
            
//             file << std::scientific << std::setprecision(8);
//             file << current_time << "\t"
//                 << this->accum_wall_loss[0] * this->q_times_wp/time_diff << "\t"
//                 << this->accum_wall_loss[1] * this->q_times_wp/time_diff << "\t"
//                 << this->accum_energy_loss[0] * this->mass * this->weight * 0.5 / time_diff << "\t"
//                 << this->accum_energy_loss[1] * this->mass * this->weight * 0.5 / time_diff << "\t"
//                 << this->total_number_particles << "\t"
//                 << KE_ave
//                 <<"\n";

//             file.close();
//             std::cout << "Number of " << this->name << " is " << total_number_particles << std::endl;
//             std::cout << "  " << std::endl;
//         } else {
//             std::ofstream file(dir_name + "/ParticleAveDiagnostic_" + this->name + ".dat");
//             if (!file) {
//                 std::cerr << "Error opening file\n";
//                 return;
//             }
//             file << "Left curr (A/m^2), right curr (A/m^2), left power (W/m^2), right power (W/m^2) \n";
//             file << std::scientific << std::setprecision(8);
//             file << this->accum_wall_loss[0] * this->q_times_wp/time_diff << "\t"
//                 << this->accum_wall_loss[1] * this->q_times_wp/time_diff << "\t"
//                 << this->accum_energy_loss[0] * this->mass * this->weight * 0.5 / time_diff << "\t"
//                 << this->accum_energy_loss[1] * this->mass * this->weight * 0.5 / time_diff 
//                 <<"\n";

//             file.close();
//         }
//     }

// }



std::vector<charged_particle> read_charged_particle_inputs(const std::string& directory_path, const domain& world){
    
    std::vector<charged_particle> particle_list;
    std::vector<double> mass_in, charge_in, n_ave, temp_in, v_drift;
    std::vector<size_t> num_part_thread, factor, EEDF_type;
    std::vector<std::string> particle_names;
    std::vector<int> index_order, number_space_coordinates, number_velocity_coordinates;
    int count_number_particles = 0;
    if (mpi_vars::mpi_rank == 0) {
        std::cout << " "<< std::endl;
        std::cout << "Reading charged particle inputs "<< std::endl;
        std::cout << "---------------------------------------- "<< std::endl;
    }
    for (int rank_num = 0; rank_num < mpi_vars::mpi_size; rank_num++){
        if (mpi_vars::mpi_rank == rank_num) {
            // Open the directory
            DIR* dir = opendir(directory_path.c_str());
            if (!dir) {
                perror("opendir");
                exit(EXIT_FAILURE);
            }

            struct dirent* entry;
            while ((entry = readdir(dir)) != nullptr) {
                if (entry->d_type == DT_REG) {  // regular file
                    count_number_particles++;
                    std::string filename = directory_path + entry->d_name;
                    std::string line;
                    std::ifstream file(filename);
                    if (!file) {
                        std::cerr << "Error: Unable to open file " << filename << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    std::getline(file, line);
                    std::istringstream iss(line);
                    std::string name;
                    iss >> name;
                    particle_names.push_back(name);
                    iss.clear();
                    std::getline(file, line);
                    iss.str(line);
                    double charge;
                    iss >> charge;
                    charge_in.push_back(charge * constants::elementary_charge);
                    iss.clear();
                    std::getline(file, line);
                    iss.str(line);
                    double mass;
                    iss >> mass;
                    mass = mass * constants::mass_amu;
                    if (std::abs(constants::electron_mass - mass ) /constants::electron_mass < 1e-3 ) {
                        mass = constants::electron_mass;
                    } else {
                        mass = mass - charge * constants::electron_mass; // assume put in neutral mass, so subtract electron mass for momentum/energy conservation in collisions
                    }
                    mass_in.push_back(mass);
                    iss.clear();
                    std::getline(file, line);
                    iss.str(line);
                    size_t number_part;
                    iss >> number_part;
                    num_part_thread.push_back(number_part);
                    iss.clear();
                    std::getline(file, line);
                    iss.str(line);
                    double n_ave_temp;
                    iss >> n_ave_temp;
                    n_ave.push_back(n_ave_temp);
                    iss.clear();
                    std::getline(file, line);
                    iss.str(line);
                    double temp;
                    iss >> temp;
                    temp_in.push_back(temp);
                    iss.clear();
                    std::getline(file, line);
                    iss.str(line);
                    double final_fac;
                    iss >> final_fac;
                    factor.push_back(size_t(final_fac * number_part));
                    iss.clear();
                    std::getline(file, line);
                    iss.str(line);
                    int EEDF_type_temp;
                    iss >> EEDF_type_temp;
                    EEDF_type.push_back(EEDF_type_temp);
                    iss.clear();
                    std::getline(file, line);
                    iss.str(line);
                    double drift;
                    iss >> drift;
                    v_drift.push_back(drift);
                    iss.clear();
                    std::getline(file, line);
                    iss.str(line);
                    int number_space, number_velocity;
                    iss >> number_space >> number_velocity;
                    number_space_coordinates.push_back(number_space);
                    number_velocity_coordinates.push_back(number_velocity);
                    iss.clear();
                    file.close();
                }
            }   
        
            closedir(dir);

            // Initialize index_order with indices [0, 1, 2, ..., count_number_particles - 1]
            index_order.resize(count_number_particles);
            std::iota(index_order.begin(), index_order.end(), 0);

            
            // Sort indices based on charge-to-mass ratio (q/m)
            std::sort(index_order.begin(), index_order.end(), [&](int i, int j) {
                return (charge_in[i] / mass_in[i]) < (charge_in[j] / mass_in[j]);
            });
        
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    

    //Create the particles in the sorted order 
    for (int num_part = 0; num_part < count_number_particles; num_part++){
        int i = index_order[num_part];
        std::string name = particle_names[i];
        double mass = mass_in[i];
        double charge = charge_in[i];
        size_t number_in = num_part_thread[i];
        size_t final_in = factor[i];
        int number_space = number_space_coordinates[i];
        int number_velocity = number_velocity_coordinates[i];
        charged_particle temp_particle(mass, charge, number_in, final_in, name, world.number_nodes);
        
        temp_particle.initialize_number_coordinates(number_space, number_velocity);
        temp_particle.initialize_weight(n_ave[i], world.length_domain);
        temp_particle.initialize_rand_maxwellian(temp_in[i], v_drift[i]);
        temp_particle.initialize_rand_position_uniform(world);
        particle_list.push_back(temp_particle);
    }

    for (int i = 0; i < particle_list.size(); i++) {
        particle_list[i].print_out();
    }

    return particle_list;
    

}


