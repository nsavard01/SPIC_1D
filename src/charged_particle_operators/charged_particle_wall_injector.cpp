
#include "charged_particle_operators/charged_particle_wall_injector.hpp"
#include "rand_gen/maxwell_generator.hpp"


charged_particle_wall_injector::charged_particle_wall_injector(std::vector<charged_particle>& particle_list, std::vector<int>& particle_indx, std::vector<double>& current_density, std::vector<std::vector<double>>& v_3D, 
    std::vector<double>& v_therm, std::vector<double>& particle_location) {
    this->particle_indx = particle_indx;
    this->current_density = current_density;
    this->v_therm = v_therm; 
    this->v_3D = v_3D;
    int total_thread_count = omp_get_max_threads() * mpi_vars::mpi_size;
    this->accumulated_number_particles.resize(this->particle_indx.size());
    this->particle_location = particle_location;
    for (int i = 0; i < this->current_density.size(); i++) {
        // Divide current density by number of threads, so have equal amount of particle introduced in each thread
        this->current_density[i] = this->current_density[i]/double(total_thread_count);
        this->accumulated_number_particles[i].resize(omp_get_max_threads());
        for (int i_thread = 0; i_thread < omp_get_max_threads(); i_thread++) {
            this->accumulated_number_particles[i][i_thread].resize(1, 0);
        }
    }

    for (int inj_indx = 0; inj_indx < this->particle_indx.size(); inj_indx++) {
        int part_indx = this->particle_indx[inj_indx];
        charged_particle& particle = particle_list[part_indx];
        if (particle.number_particles_injected.empty()) {
            particle.number_particles_injected.resize(omp_get_max_threads());
            particle.time_step_injected.resize(omp_get_max_threads());
            particle.number_unique_injections = 0;
            for (int i_thread = 0; i_thread < omp_get_max_threads(); i_thread++) {
                particle.time_step_injected[i_thread].push_back({});
            }
        }
        for (int i_thread = 0; i_thread < omp_get_max_threads(); i_thread++) {
            particle.number_particles_injected[i_thread].push_back(0);
            particle.time_step_injected[i_thread][particle.number_unique_injections].resize(100, 0.0); // set initial number guess twice the expected amount , in case ever del_t increase
        }
        particle.number_unique_injections++;
    }

}

void charged_particle_wall_injector::setup_diagnostics(const std::string& dir_name, std::vector<charged_particle>& particle_list) {
    
    std::vector<int> amount_injections(particle_list.size(), 0);

    for (int inj_indx = 0; inj_indx < this->particle_indx.size(); inj_indx++) {
        int part_indx = this->particle_indx[inj_indx];
        charged_particle& particle = particle_list[part_indx];
        int number_injection = amount_injections[part_indx];
        std::ofstream file(dir_name + "/charged_particles/" + particle.name + "/operations/wall_injection_properties_" + std::to_string(number_injection) + ".dat");
        double total_current_density = this->current_density[inj_indx] * omp_get_max_threads() * mpi_vars::mpi_size;
        double Temp = this->v_therm[inj_indx] * this->v_therm[inj_indx] / std::abs(particle.q_over_m);
        // Write header (optional)
        file << "Current density (A/m^2), thermal temperature (eV), v_x (m/s), v_y (m/s), v_z (m/s), wall node \n";
        file << std::scientific << std::setprecision(8);
        file << total_current_density << "\t"
            << Temp << "\t"
            << this->v_3D[inj_indx][0] << "\t"
            << this->v_3D[inj_indx][1] << "\t"
            << this->v_3D[inj_indx][2]<< "\t"
            << std::round(this->particle_location[inj_indx]) << "\n";


        file.close();

        file.open(dir_name + "/charged_particles/" + particle.name + "/operations/wall_injection_diagnostics_" + std::to_string(number_injection) + ".dat");
        file << "Number particles released \n";

        file.close();

        amount_injections[part_indx]++;
    }
}

void charged_particle_wall_injector::write_average_diagnostics(const std::string& dir_name, const std::vector<charged_particle>& particle_list) {
    
    std::vector<size_t> total_amount_injected(this->particle_indx.size(),0);
    for (int inj_indx = 0; inj_indx < this->particle_indx.size(); inj_indx++) {
        for (int i_thread = 0; i_thread < omp_get_max_threads(); i_thread++) {
            total_amount_injected[inj_indx] += this->accumulated_number_particles[inj_indx][i_thread][0];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, total_amount_injected.data(), this->particle_indx.size(), mpi_vars::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    
    if (mpi_vars::mpi_rank == 0) {
        std::vector<int> amount_injections(particle_list.size(), 0);

        for (int inj_indx = 0; inj_indx < this->particle_indx.size(); inj_indx++) {
            int part_indx = this->particle_indx[inj_indx];
            const charged_particle& particle = particle_list[part_indx];
            int number_injection = amount_injections[part_indx];
            std::ofstream file(dir_name + "/charged_particles/" + particle.name + "/operations/wall_injection_diagnostics_average_" + std::to_string(number_injection) + ".dat");
            file << "Number particles released \n";
            file << total_amount_injected[inj_indx] << "\n";


            file.close();
        }
    }
}

void charged_particle_wall_injector::write_diagnostics(const std::string& dir_name, const std::vector<charged_particle>& particle_list) {
    
    std::vector<size_t> total_amount_injected(this->particle_indx.size(),0);
    for (int inj_indx = 0; inj_indx < this->particle_indx.size(); inj_indx++) {
        for (int i_thread = 0; i_thread < omp_get_max_threads(); i_thread++) {
            total_amount_injected[inj_indx] += this->accumulated_number_particles[inj_indx][i_thread][0];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, total_amount_injected.data(), this->particle_indx.size(), mpi_vars::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    
    if (mpi_vars::mpi_rank == 0) {
        std::vector<int> amount_injections(particle_list.size(), 0);

        for (int inj_indx = 0; inj_indx < this->particle_indx.size(); inj_indx++) {
            int part_indx = this->particle_indx[inj_indx];
            const charged_particle& particle = particle_list[part_indx];
            int number_injection = amount_injections[part_indx];
            std::ofstream file(dir_name + "/charged_particles/" + particle.name + "/operations/wall_injection_diagnostics_" + std::to_string(number_injection) + ".dat", std::ios::app);
            file << std::scientific << std::setprecision(8);
            file << total_amount_injected[inj_indx] << "\n";


            file.close();
        }
    }
}

void charged_particle_wall_injector::reset_diagnostics() {
    for (int i = 0; i < this->current_density.size(); i++) {
        for (int i_thread = 0; i_thread < omp_get_max_threads(); i_thread++) {
            this->accumulated_number_particles[i][i_thread][0] = 0;
        }
    }
    
}

void charged_particle_wall_injector::print_out() {
    if (mpi_vars::mpi_rank == 0) {
        std::cout  << std::endl;
        std::cout << " ------------ " << std::endl;
        std::cout << "Particle Injection Operator " << std::endl;
        for (int i = 0; i < this->particle_indx.size(); i++) {
            std::cout << "Particle index " << this->particle_indx[i] << std::endl;
            std::cout << "With current density: " << this->current_density[i] << " A/m^2 per thread" << std::endl;
            std::cout << "Particle location: " << this->particle_location[i] << std::endl;
            std::cout << "thermal velocity: " << this->v_therm[i] << std::endl;
            std::cout << "v_x: " << v_3D[i][0] <<  ", v_y: " << v_3D[i][1] << " , v_z: " << v_3D[i][2] << std::endl;
            std::cout << std::endl;
        }
        std::cout << " ------------ " << std::endl;
        std::cout  << std::endl;
    }
}


void charged_particle_wall_injector::run(const int thread_id, const double current_time, const double del_t, std::vector<charged_particle>& particle_list, const domain& world) {
    std::vector<int> amount_injections(particle_list.size(), 0);
    for (int i = 0; i < this->particle_indx.size(); i++){
        int part_indx = this->particle_indx[i];
        charged_particle& part = particle_list[part_indx];
        const double particle_position = this->particle_location[i];
        const double v_x_drift = this->v_3D[i][0];
        const double v_y_drift = this->v_3D[i][1];
        const double v_z_drift = this->v_3D[i][2];
        const double v_therm_local = this->v_therm[i] * ((v_x_drift>0) - (v_x_drift<0));
        double number_particles_inject = std::abs(del_t * this->current_density[i]/part.q_times_wp); 
        size_t number_selected = size_t(number_particles_inject); 
        if (pcg32_random_r() < (number_particles_inject - number_selected)) {
            number_selected++;
        }
        int inject_indx = amount_injections[part_indx];
        amount_injections[part_indx]++;
        part.number_particles_injected[thread_id][inject_indx] = number_selected;
        if (part.time_step_injected[thread_id][inject_indx].size() < number_selected) {
            part.time_step_injected[thread_id][inject_indx].resize(number_selected + 1, 0);
        }
        size_t number_particles = part.number_particles[thread_id][0];
        // references to particle arrays
        std::vector<double>& xi_local = part.xi[thread_id];
        std::vector<double>& v_x_local = part.v_x[thread_id];
        bool use_vy = (part.number_velocity_coordinates > 1);
        bool use_vz = (part.number_velocity_coordinates > 2);
        double v_x_temp, v_y_temp, v_z_temp;
        for (size_t part_indx = 0; part_indx < number_selected; part_indx++){
            maxwellian_3D_flux(v_x_temp, v_y_temp, v_z_temp, v_therm_local, v_x_drift);
            xi_local[number_particles] = particle_position;
            v_x_local[number_particles] = v_x_temp;
            part.time_step_injected[thread_id][inject_indx][part_indx] = pcg32_random_r() * del_t;
            if (use_vy) {
                part.v_y[thread_id][number_particles] = v_y_temp + v_y_drift;
            }
            if (use_vz) {
                part.v_z[thread_id][number_particles] = v_z_temp + v_z_drift;
            }
            number_particles++;
        }
        
        this->accumulated_number_particles[i][thread_id][0] += number_selected;
        
    }
    

}


