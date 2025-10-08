#include "collisions/null_collider.hpp"
#include "globals/constants.hpp"
#include "globals/mpi_vars.hpp"
#include "rand_gen/maxwell_generator.hpp"
#include <regex>
#include <cstdlib>
#include <omp.h>
#include <algorithm>
#include <numeric>
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <dirent.h>
#include <fstream>
#include <sstream>


static double sin_third_rot = std::sin(2.0 * M_PI / 3.0);
static double cos_third_rot = std::cos(2.0 * M_PI / 3.0);

null_collider::null_collider(int primary_idx, int number_targets, const std::vector<int>& target_idx, const std::vector<int>& number_collisions_per_target,  
    const std::vector<std::vector<std::vector<double>>> &sigma_array, 
    const std::vector<double> &energy_array, const std::vector<std::vector<double>> &energy_threshold,
    const std::vector<std::vector<int>> &collision_type_per_target, const std::vector<std::vector<std::vector<int>>> &product_indices,
    const std::vector<double>& reduced_mass, const std::vector<double>& reduced_mass_ionization)
    {
    this->number_collisions_per_target = number_collisions_per_target;
    this->primary_idx = primary_idx;
    this->number_targets = number_targets;
    this->target_idx = target_idx;
    this->sigma_array = sigma_array;
    this->energy_array = energy_array;
    this->energy_threshold = energy_threshold;
    this->collision_type_per_target = collision_type_per_target;
    this->product_indices = product_indices;
    this->reduced_mass = reduced_mass;
    this->reduced_mass_ionization = reduced_mass_ionization;
    this->total_amount_collidable_particles = 0;
    if (number_targets > 0) {
        this->total_incident_energy.resize(this->number_targets);
        this->total_energy_loss.resize(this->number_targets);
        this->total_amount_collisions.resize(this->number_targets);
        this->collision_id.resize(this->number_targets);
        for (int i = 0; i < this->number_targets; i++) {
            int num_collisions = this->number_collisions_per_target[i];
            this->total_incident_energy[i].resize(num_collisions, 0.0);
            this->total_energy_loss[i].resize(num_collisions, 0.0);
            this->total_amount_collisions[i].resize(num_collisions, 0);
            this->collision_id[i].resize(num_collisions);
            for (int j = 0; j < this->number_collisions_per_target[i]; j++) {
                this->collision_id[i][j] = j;
            }
        }
        int max_threads = omp_get_max_threads();
        this->total_incident_energy_thread.resize(max_threads);
        this->total_amount_collisions_thread.resize(max_threads);
        this->total_energy_loss_thread.resize(max_threads);
        for (int i_thread = 0; i_thread < max_threads; i_thread++) {
            this->total_incident_energy_thread[i_thread].resize(this->number_targets);
            this->total_amount_collisions_thread[i_thread].resize(this->number_targets);
            this->total_energy_loss_thread[i_thread].resize(this->number_targets);
            for (int i = 0; i < this->number_targets; i++) {
                int num_collisions = this->number_collisions_per_target[i];
                this->total_incident_energy_thread[i_thread][i].resize(num_collisions, 0);
                this->total_amount_collisions_thread[i_thread][i].resize(num_collisions, 0);
                this->total_energy_loss_thread[i_thread][i].resize(num_collisions, 0);
            }
        }
    }
    
    
 
}

void null_collider::set_null_frequency(const std::vector<charged_particle>& particle_list, const std::vector<target_particle>& target_particle_list) {
    if (this->number_targets > 0) {
        double primary_mass = particle_list[this->primary_idx].mass;
        std::vector<double> target_density(this->number_targets);
        std::vector<double> red_mass(this->number_targets);
        for (int i = 0; i < this->number_targets; i++){
            target_density[i] = target_particle_list[this->target_idx[i]].average_density;
            double target_mass = target_particle_list[this->target_idx[i]].mass;
            red_mass[i] = (primary_mass * target_mass)/(primary_mass + target_mass);
        }
        this->null_frequency = 0.0;
        double v_r;
        int size_array = this->energy_array.size();
        for (int i = 0; i < size_array; i++) {
            double freq_sum = 0.0;
            double E = this->energy_array[i];
            for (int t_idx = 0; t_idx < this->number_targets; t_idx++) {
                v_r = std::sqrt(2.0 * E * constants::elementary_charge/red_mass[t_idx]);
                for (int coll_idx = 0; coll_idx < this->number_collisions_per_target[t_idx]; coll_idx++){
                    freq_sum += this->sigma_array[t_idx][coll_idx][i] * target_density[t_idx] * v_r;
                }
            }
            
            this->null_frequency = std::max(this->null_frequency, freq_sum);
        }

    }
}

void null_collider::print_out(const std::vector<charged_particle>& particle_list, const std::vector<target_particle>& target_particle_list) const {
    if (mpi_vars::mpi_rank == 0 && this->number_targets > 0) {
        std::cout << "------------------------------ " << std::endl;
        std::cout << "Primary particle is " << particle_list[this->primary_idx].name << std::endl;
        std::cout << "Amount of targets is " << this->number_targets << std::endl;
        std::cout << "Null frequency: " << this->null_frequency << std::endl;
        std::cout << "----- " << std::endl;
        std::cout << std::endl;
        for (int l=0; l<this->number_targets;l++){
            std::cout << "Target # " << l << ": " << target_particle_list[this->target_idx[l]].name << std::endl;
            std::cout << "Reduced mass " << this->reduced_mass[l] << " kg" << std::endl;
            std::cout << "Reduced mass ionization " << this->reduced_mass_ionization[l] << " kg" << std::endl;
            for (int coll_idx = 0; coll_idx < this->number_collisions_per_target[l]; coll_idx++){
                std::cout << "\t" << "Collision # " << coll_idx << std::endl;
                std::cout << "\t" << "Threshold energy (eV) " << this->energy_threshold[l][coll_idx] << std::endl;
                int coll_type = this->collision_type_per_target[l][coll_idx];
                std::cout << "\t" << "Collision type is ";
                if (coll_type == 1) {
                    std::cout << "elastic" << std::endl;
                    std::cout << "\t" << "Products are: " << particle_list[this->product_indices[l][coll_idx][0]].name << " " << target_particle_list[this->product_indices[l][coll_idx][1]].name << std::endl;
                } else if (coll_type == 2) {
                    std::cout << "ionization" << std::endl;
                    std::cout << "\t" << "Products are: " << particle_list[this->product_indices[l][coll_idx][0]].name << " " << particle_list[this->product_indices[l][coll_idx][1]].name << "  e " << std::endl;
                } else if (coll_type == 3) {
                    std::cout << "excitation" << std::endl;
                    std::cout << "\t" << "Products are: " << particle_list[this->product_indices[l][coll_idx][0]].name << " " << target_particle_list[this->product_indices[l][coll_idx][1]].name << std::endl;
                } else if (coll_type == 4) {
                    std::cout << "charge exchange" << std::endl;
                    std::cout << "\t" << "Products are: " << particle_list[this->product_indices[l][coll_idx][0]].name << " " << target_particle_list[this->product_indices[l][coll_idx][1]].name << std::endl;
                }
            
                std::cout << std::endl;
            }
        }
        std::cout << "------------------------------ " << std::endl;
    }
}

inline void double_product_isotropic(const double &primary_mass, const double &target_mass, const double &reduced_mass, const double &del_E, double (&incident_velocity)[3], double (&target_velocity)[3]) {
    
    double e_vector[3];
    double cos_theta, speed_per_particle, phi, sin_theta, cos_phi, sin_phi, P_beginning, V_cm;
    int i;
    double mass_sum = primary_mass + target_mass;

    speed_per_particle = std::sqrt(2.0 * del_E * reduced_mass / primary_mass/primary_mass);

    cos_theta = 1.0 - 2.0 * pcg32_random_r();
    phi = pcg32_random_r() * 2.0 * M_PI;
    sin_theta = std::sqrt(1.0 - cos_theta*cos_theta);
    cos_phi = std::cos(phi);
    sin_phi = std::sin(phi);
    

    e_vector[0] = cos_phi * sin_theta;
    e_vector[1] = sin_phi * sin_theta;
    e_vector[2] = cos_theta;
    
    for (i=0;i<3;i++){
        P_beginning = primary_mass * incident_velocity[i] + target_mass * target_velocity[i];
        V_cm = P_beginning/mass_sum;
        incident_velocity[i] = e_vector[i] * speed_per_particle + V_cm;
        target_velocity[i] = (P_beginning - primary_mass * incident_velocity[i])/target_mass;
    }

}

inline void triple_product_isotropic(const double &primary_mass, const double &ion_mass, const double &target_mass, const double &reduced_mass, const double &del_E, 
    double (&incident_velocity)[3], double (&target_velocity)[3], double (&third_velocity)[3]) {
    
    double e_vector[3], y_vector[3], u_vector[3];
    double cos_theta, speed_per_particle, phi, sin_theta, cos_phi, sin_phi, cos_theta_new, sin_theta_new, speed_electron, P_beginning, V_cm;
    int i;
    double mass_sum = primary_mass + target_mass;


    speed_per_particle = std::sqrt(2.0 * del_E * reduced_mass / primary_mass/primary_mass);

    cos_theta = 1.0 - 2.0 * pcg32_random_r();
    phi = pcg32_random_r() * 2.0 * M_PI;
    sin_theta = std::sqrt(1.0 - cos_theta*cos_theta);
    cos_phi = std::cos(phi);
    sin_phi = std::sin(phi);
    e_vector[0] = cos_phi * sin_theta;
    e_vector[1] = sin_phi * sin_theta;
    e_vector[2] = cos_theta;

    cos_theta_new = cos_third_rot * cos_theta - sin_third_rot * sin_theta;
    sin_theta_new = cos_theta * sin_third_rot + cos_third_rot * sin_theta;
    y_vector[0] = cos_phi * sin_theta_new;
    y_vector[1] = sin_phi * sin_theta_new;
    y_vector[2] = cos_theta_new;

    phi = pcg32_random_r() * 2.0 * M_PI;
    cos_phi = std::cos(phi);
    sin_phi = std::sin(phi);

    u_vector[0] = e_vector[1] * y_vector[2] - e_vector[2] * y_vector[1];
    u_vector[1] = e_vector[2] * y_vector[0] - e_vector[0] * y_vector[2];
    u_vector[2] = e_vector[0] * y_vector[1] - e_vector[1] * y_vector[0];

    speed_electron = std::sqrt(2.0 * del_E * reduced_mass/ constants::electron_mass/constants::electron_mass);

    for (i=0;i<3;i++){
        u_vector[i] = -e_vector[i] * cos_third_rot * (cos_phi - 1.0) + y_vector[i] * cos_phi + u_vector[i] * sin_phi;
        P_beginning = primary_mass * incident_velocity[i] + target_mass * target_velocity[i];
        V_cm = P_beginning/mass_sum;
        incident_velocity[i] = e_vector[i] * speed_per_particle + V_cm;
        third_velocity[i] = u_vector[i] * speed_electron + V_cm;
        target_velocity[i] = (P_beginning - primary_mass * incident_velocity[i] - constants::electron_mass * third_velocity[i])/ion_mass;
    }


}

// void Null_Collision::initialize_data_files(const std::string& dir_name, std::vector<Particle>& particle_list, std::vector<Target_Particle>& target_particle_list) const {

//     if (Constants::mpi_rank == 0) {
//         std::string binary_folder = dir_name + "/BinaryCollisions/" + particle_list[this->primary_idx].name + "_on_" + target_particle_list[this->target_idx].name;
//         if (!createDirectory(binary_folder)) {
//             std::cerr << "Save directory not successfully created!" << std::endl;
//             MPI_Abort(MPI_COMM_WORLD, 1);
//         }
//         std::ofstream file(binary_folder + "/CollisionProperties.dat");
//         if (!file) {
//             std::cerr << "Error opening file \n";
//             return;
//         }

//         file << "Coll #, collType, E_thres (eV), maxSigma (m^2), EatMaxSigma (eV) \n";
//         file << std::scientific << std::setprecision(8);
        
//         for (int i=0; i< this->number_collisions;i++){
//             size_t max_indx;
//             double max_val = 0.0;
//             for (int j= 0; j< this->length_arrays;j++){
//                 if (this->sigma_array[i][j] > max_val) {
//                     max_indx = j;
//                     max_val = this->sigma_array[i][j];
//                 }
//             }
//             file << i << "\t"
//                 << this->collision_type[i] << "\t"
//                 << this->energy_threshold[i] << "\t"
//                 << this->sigma_array[i][max_indx] << "\t"
//                 << this->energy_array[max_indx]
//                 <<"\n";
//         }
        
//         file.close();
//         for (int i=0; i< this->number_collisions;i++){
//             file.open(binary_folder + "/CollisionDiag_" + std::to_string(i+1) + ".dat");
//             file << "CollRatio, AveEnergyLoss (eV), AveIncidentEnergy (eV), P_loss(W/m^2), aveCollFreq (Hz/m^2) \n";
//             file.close();
//         }
        
//     }
// }


// void Null_Collision::gather_mpi() {
//     MPI_Allreduce(MPI_IN_PLACE, this->total_incident_energy.data(), this->number_collisions, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//     MPI_Allreduce(MPI_IN_PLACE, this->total_energy_loss.data(), this->number_collisions, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//     MPI_Allreduce(MPI_IN_PLACE, this->total_amount_collisions.data(), this->number_collisions, Constants::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
//     MPI_Allreduce(MPI_IN_PLACE, &this->total_amount_collidable_particles, 1, Constants::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
// }

// void Null_Collision::diag_write(const std::string& dir_name, std::vector<Particle>& particle_list, std::vector<Target_Particle>& target_particle_list, const double& time_diff, bool average_bool) {
//     if (Constants::mpi_rank == 0) {
//         std::string binary_folder = dir_name + "/BinaryCollisions/" + particle_list[this->primary_idx].name + "_on_" + target_particle_list[this->target_idx].name;
//         if (!average_bool) {   
//             for (int i=0; i< this->number_collisions;i++){
//                 std::ofstream file(binary_folder + "/CollisionDiag_" + std::to_string(i+1) + ".dat", std::ios::app);
//                 file << std::scientific << std::setprecision(8);
//                 file << double(this->total_amount_collisions[i])/double(this->total_amount_collidable_particles) << "\t"
//                     << this->total_energy_loss[i] * 0.5 / Constants::elementary_charge/double(this->total_amount_collisions[i]) << "\t"
//                     << this->total_incident_energy[i] * particle_list[this->primary_idx].mass * 0.5 / Constants::elementary_charge/double(this->total_amount_collisions[i]) << "\t"
//                     << this->total_energy_loss[i] * 0.5 * particle_list[this->primary_idx].weight / time_diff << "\t"
//                     << double(this->total_amount_collisions[i]) * particle_list[this->primary_idx].weight / time_diff
//                     << "\n";
//                 file.close();
//             }
//         } else {
//             std::ofstream file(binary_folder + "/AveCollisionDiag.dat");
//             file << "Coll #, CollRatio, AveEnergyLoss (eV), AveIncidentEnergy (eV), P_loss(W/m^2), aveCollFreq (Hz/m^2) \n";
//             file << std::scientific << std::setprecision(8);
//             for (int i=0; i< this->number_collisions;i++){
//                 file << i+1 << "\t"
//                     << double(this->total_amount_collisions[i])/double(this->total_amount_collidable_particles) << "\t"
//                     << this->total_energy_loss[i] * 0.5 / Constants::elementary_charge/double(this->total_amount_collisions[i]) << "\t"
//                     << this->total_incident_energy[i] * particle_list[this->primary_idx].mass * 0.5 / Constants::elementary_charge/double(this->total_amount_collisions[i]) << "\t"
//                     << this->total_energy_loss[i] * 0.5 * particle_list[this->primary_idx].weight / time_diff << "\t"
//                     << double(this->total_amount_collisions[i]) * particle_list[this->primary_idx].weight / time_diff
//                     << "\n";
//             }
//             file.close();
//         }
        
//     }
// }



void null_collider::generate_null_collisions(const int thread_id, std::vector<charged_particle> &particle_list, const std::vector<target_particle> &target_particle_list, const double time_step){
    
    if (this->number_targets > 0) {
        
        // initialize local variables
        const int length_energy_array = this->energy_array.size();
        const int number_targets_local = this->number_targets;
        const double min_energy = this->energy_array[0];
        const double max_energy = this->energy_array.back();
        const double null_frequency_local = this->null_frequency;
        const std::vector<int>& target_idx_local = this->target_idx;
        const std::vector<double>& energy_array_local = this->energy_array;
        const std::vector<std::vector<std::vector<double>>>& sigma_array_local = this->sigma_array;
        const std::vector<std::vector<double>>& energy_threshold_local = this->energy_threshold;
        const std::vector<double>& reduced_mass_local = this->reduced_mass;
        const std::vector<double>& reduced_mass_ionization_local = this->reduced_mass_ionization;
        const std::vector<int>& number_collisions_per_target_local = this->number_collisions_per_target;
        const std::vector<std::vector<std::vector<int>>>& product_indices_local = this->product_indices;
        const std::vector<std::vector<int>>& collision_type_per_target_local = this->collision_type_per_target;
        std::vector<std::vector<size_t>>& total_collisions = this->total_amount_collisions_thread[thread_id];
        std::vector<std::vector<double>>& energy_loss = this->total_energy_loss_thread[thread_id];
        std::vector<std::vector<double>>& tot_incident_energy = this->total_incident_energy_thread[thread_id];
        std::vector<double> target_mass(number_targets_local),

        // generate local vectors for target and diagnostics
        target_density(number_targets_local), v_therm(number_targets_local);
        for (int i = 0; i < number_targets_local; i++) {
            int idx = target_idx_local[i];
            target_mass[i] = target_particle_list[idx].mass;
            target_density[i] = target_particle_list[idx].average_density;
            v_therm[i] = target_particle_list[idx].v_therm;
        }

        // Calculate 
        double P_null = 1.0 - std::exp(-null_frequency_local * time_step);
        if (P_null > 0.05){
            std::cout << "P_null greater than 5% " << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        charged_particle& primary_particle = particle_list[primary_idx];
        size_t particle_indx, end_indx;
        double particle_location, incident_velocity[3], target_velocity[3], velocity_CM[3], speed_CM, energy_CM, interp_d,
            freq_sum, incident_energy, del_E, test_freq;
        const size_t initial_amount_collidable_particles = primary_particle.number_collidable_particles[thread_id][0];
        const double primary_mass = primary_particle.mass;
        size_t number_total_particles = initial_amount_collidable_particles;
        double number_selected_real = P_null * static_cast<double>(number_total_particles);
        size_t number_selected = static_cast<size_t>(number_selected_real);
        double Rand = pcg32_random_r();
        if (Rand < (number_selected_real - number_selected)) {
            number_selected++;
        }
        std::vector<double>& xi_local = primary_particle.xi[thread_id];
        std::vector<double>& v_x_local = primary_particle.v_x[thread_id];
        std::vector<double>& v_y_local = primary_particle.v_y[thread_id];
        std::vector<double>& v_z_local = primary_particle.v_z[thread_id];
        int indx_low = 0, indx_high = length_energy_array-1, indx_middle;
        double energy_low, energy_high;
        energy_low = energy_array_local[indx_low];
        energy_high = energy_array_local[indx_high];
        for (size_t part_num=0;part_num<number_selected;part_num++){
            particle_indx = static_cast<size_t>(number_total_particles * pcg32_random_r());
            particle_location = xi_local[particle_indx];
            incident_velocity[0] = v_x_local[particle_indx];
            incident_velocity[1] = v_y_local[particle_indx];
            incident_velocity[2] = v_z_local[particle_indx];
            bool collided = false;
            freq_sum = 0.0;
            test_freq = pcg32_random_r() * null_frequency_local; // Rand collision freq for determining 
            for (int t_idx = 0; t_idx < number_targets_local; t_idx++) {
                double target_density_local = target_density[t_idx];
                double red_mass = reduced_mass_local[t_idx];
                double target_mass_local = target_mass[t_idx];
                maxwellian_3D(target_velocity[0], target_velocity[1], target_velocity[2], v_therm[t_idx], 0.0);
                speed_CM = 0.0;
                for (int iter=0;iter<3;iter++){
                    velocity_CM[iter] = incident_velocity[iter] - target_velocity[iter];
                    speed_CM += velocity_CM[iter] * velocity_CM[iter];
                }
                energy_CM = speed_CM * 0.5 * red_mass / constants::elementary_charge; // CM energy in eV

                // Make sure its in range
                if (energy_CM <= min_energy){
                    // Use minimum sigma
                    indx_low = 0;
                    interp_d = 0.0;
                } else if (energy_CM >= max_energy) {
                    // Use maximum sigma
                    indx_high = length_energy_array-1;
                    interp_d = 1.0;
                } else {
                    if (energy_CM < energy_low) {
                        indx_low = 0;
                    } else if (energy_CM > energy_high) {
                        indx_high = length_energy_array-1;
                    }
                    while (indx_low != indx_high-1) {
                        indx_middle = (indx_low + indx_high)/2;
                        if (energy_array_local[indx_middle] < energy_CM) {
                            indx_low = indx_middle;
                        } else if (energy_array_local[indx_middle] > energy_CM) {
                            indx_high = indx_middle;
                        } else {
                            indx_low = indx_middle;
                            indx_high = indx_low + 1;
                        }
                    }
                    energy_low = energy_array_local[indx_low];
                    energy_high = energy_array_local[indx_high];
                    interp_d = (energy_CM - energy_low)/(energy_high - energy_low);
                }
                speed_CM = std::sqrt(speed_CM); // convert to physical relative speed
                int tot_num_collisions = number_collisions_per_target_local[t_idx];
                for (int coll_idx=0;coll_idx<tot_num_collisions;coll_idx++){
                    double thres_E = energy_threshold_local[t_idx][coll_idx];
                    if (energy_CM > thres_E) {
                        freq_sum += (sigma_array_local[t_idx][coll_idx][indx_low] * (1.0 - interp_d) + sigma_array_local[t_idx][coll_idx][indx_high] * interp_d) * speed_CM * target_density_local;
                        collided = test_freq <= freq_sum;
                        if (collided) {
                            // Collision selected
                            total_collisions[t_idx][coll_idx]++;
                            incident_energy = incident_velocity[0]*incident_velocity[0] + incident_velocity[1]*incident_velocity[1] + incident_velocity[2]*incident_velocity[2];
                            tot_incident_energy[t_idx][coll_idx] += incident_energy;
                            del_E = (energy_CM - thres_E) * constants::elementary_charge;
                            switch (collision_type_per_target_local[t_idx][coll_idx]) {
                                case 1: {
                                    double_product_isotropic(primary_mass, target_mass_local, red_mass, del_E, incident_velocity, target_velocity);
                                    break;
                                }
                                case 2: {
                                    int secondary_product_idx = product_indices_local[t_idx][coll_idx][1]; // secondary product ion
                                    charged_particle& secondary_particle = particle_list[secondary_product_idx];
                                    charged_particle& electron_particle = particle_list[0]; // If ionization exists, electron exists at indx 0
                                    double secondary_mass = secondary_particle.mass;
                                    triple_product_isotropic(primary_mass, secondary_mass, target_mass_local, reduced_mass_ionization_local[t_idx], del_E, 
                                        incident_velocity, target_velocity, velocity_CM);
                                    // Increase number electron and ion
                                    electron_particle.number_particles[thread_id][0]++;
                                    secondary_particle.number_particles[thread_id][0]++;
    
                                    // electron set into velocity_CM
                                    size_t electron_number = electron_particle.number_particles[thread_id][0]-1;
                                    electron_particle.xi[thread_id][electron_number] = particle_location;
                                    electron_particle.v_x[thread_id][electron_number] = velocity_CM[0];
                                    electron_particle.v_y[thread_id][electron_number] = velocity_CM[1];
                                    electron_particle.v_z[thread_id][electron_number] = velocity_CM[2];
                                    
    
                                    // ion set into target_velocity
                                    size_t secondary_number = secondary_particle.number_particles[thread_id][0]-1;
                                    secondary_particle.xi[thread_id][secondary_number] = particle_location;
                                    secondary_particle.v_x[thread_id][secondary_number] = target_velocity[0];
                                    secondary_particle.v_y[thread_id][secondary_number] = target_velocity[1];
                                    secondary_particle.v_z[thread_id][secondary_number] = target_velocity[2];
                                    energy_loss[t_idx][coll_idx] += (- constants::electron_mass * (velocity_CM[0]*velocity_CM[0] +
                                        velocity_CM[1]*velocity_CM[1] + velocity_CM[2]*velocity_CM[2]) - secondary_mass * (target_velocity[0]*target_velocity[0] + 
                                        target_velocity[1]*target_velocity[1] + target_velocity[2]*target_velocity[2])); // add gain of energy in system due to introduction of target velocity
                                    break;
                                }
                                case 3: {
                                    double_product_isotropic(primary_mass, target_mass_local, red_mass, del_E, incident_velocity, target_velocity);
                                    break;
                                }
                                case 4: {
                                    incident_velocity[0] = target_velocity[0];
                                    incident_velocity[1] = target_velocity[1];
                                    incident_velocity[2] = target_velocity[2];
                                    break;
                                }
                            }
                            // Calculate energy loss and break loop
                            energy_loss[t_idx][coll_idx] += primary_mass * (incident_energy - incident_velocity[0]*incident_velocity[0] - 
                                incident_velocity[1]*incident_velocity[1] - incident_velocity[2]*incident_velocity[2]);
                            break;

                        }
                    }
                }
                if (collided) {break;} // if previous target collided, then break target loop
            }
            // switch particles with back to ensure one collision per particle
            end_indx = number_total_particles-1;
            xi_local[particle_indx] = xi_local[end_indx];
            v_x_local[particle_indx] = v_x_local[end_indx];
            v_y_local[particle_indx] = v_y_local[end_indx];
            v_z_local[particle_indx] = v_z_local[end_indx];

            xi_local[end_indx] = particle_location;
            v_x_local[end_indx] = incident_velocity[0];
            v_y_local[end_indx] = incident_velocity[1];
            v_z_local[end_indx] = incident_velocity[2];

            number_total_particles--;
        }
        primary_particle.number_collidable_particles[thread_id][0] = number_total_particles; // less collidable particles, in case have other binary collisions as well
        // sum up collision totals
        #pragma omp critical
        {
            this->total_amount_collidable_particles += initial_amount_collidable_particles;
        }
    
    }

}


void null_collider::initialize_diagnostic_files(const std::string& dir_name, const std::vector<charged_particle>& particle_list, const std::vector<target_particle>& target_particle_list) const {
    if (mpi_vars::mpi_rank == 0) {
        for (int t_idx = 0; t_idx < this->number_targets; t_idx++) {
            const charged_particle& primary_particle = particle_list[this->primary_idx];
            const target_particle& secondary_particle = target_particle_list[this->target_idx[t_idx]];
            for (int coll_idx=0; coll_idx< this->number_collisions_per_target[t_idx];coll_idx++){
                std::ofstream file(dir_name + "/charged_particles/" + primary_particle.name + "/null_collision/" + secondary_particle.name + "/collision_properties_" 
                    + std::to_string(this->collision_id[t_idx][coll_idx]) + ".dat");
                if (!file) {
                    std::cerr << "Error opening file \n";
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }

                file << "Coll #, collision type, threshold energy (eV), max sigma (m^2), energy at peak (eV) \n";
                file << std::scientific << std::setprecision(8);
                
                
                size_t max_indx;
                double max_val = 0.0;
                int length_array = this->energy_array.size();
                for (int j= 0; j< length_array;j++){
                    if (this->sigma_array[t_idx][coll_idx][j] > max_val) {
                        max_indx = j;
                        max_val = this->sigma_array[t_idx][coll_idx][j];
                    }
                }
                file << this->collision_id[t_idx][coll_idx] << "\t"
                    << this->collision_type_per_target[t_idx][coll_idx] << "\t"
                    << this->energy_threshold[t_idx][coll_idx] << "\t"
                    << this->sigma_array[t_idx][coll_idx][max_indx] << "\t"
                    << this->energy_array[max_indx]
                    <<"\n";
                file.close();

                file.open(dir_name + "/charged_particles/" + primary_particle.name + "/null_collision/" + secondary_particle.name + "/collision_diagnostics_" 
                    + std::to_string(this->collision_id[t_idx][coll_idx]) + ".dat");
                file << "Amount collisions, amount collidable,  Energy loss (J), Incident Energy (J) \n";
                file.close();
            }     
        }
    }
}

void null_collider::write_diagnostics(const std::string& dir_name, const std::vector<charged_particle>& particle_list, const std::vector<target_particle>& target_particle_list) const {
    if (mpi_vars::mpi_rank == 0) {
        for (int t_idx = 0; t_idx < this->number_targets; t_idx++) {
            const charged_particle& primary_particle = particle_list[this->primary_idx];
            const target_particle& secondary_particle = target_particle_list[this->target_idx[t_idx]];
            for (int coll_idx=0; coll_idx< this->number_collisions_per_target[t_idx];coll_idx++){
            
                std::ofstream file(dir_name + "/charged_particles/" + primary_particle.name + "/null_collision/" + secondary_particle.name + "/collision_diagnostics_" 
                    + std::to_string(this->collision_id[t_idx][coll_idx]) + ".dat", std::ios::app);
                if (!file) {
                    std::cerr << "Error opening file \n";
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                file << std::scientific << std::setprecision(8);
                file << this->total_amount_collisions[t_idx][coll_idx] << "\t" << this->total_amount_collidable_particles << "\t"
                << this->total_energy_loss[t_idx][coll_idx] * 0.5 << "\t"
                << this->total_incident_energy[t_idx][coll_idx] * 0.5 * primary_particle.mass << "\n";
                file.close();
            }     
        }
    }
}


void null_collider::write_diagnostics_average(const std::string& dir_name, const std::vector<charged_particle>& particle_list, const std::vector<target_particle>& target_particle_list) const {
    if (mpi_vars::mpi_rank == 0) {
        for (int t_idx = 0; t_idx < this->number_targets; t_idx++) {
            const charged_particle& primary_particle = particle_list[this->primary_idx];
            const target_particle& secondary_particle = target_particle_list[this->target_idx[t_idx]];
            for (int coll_idx=0; coll_idx< this->number_collisions_per_target[t_idx];coll_idx++){
            
                std::ofstream file(dir_name + "/charged_particles/" + primary_particle.name + "/null_collision/" + secondary_particle.name + "/collision_diagnostics_average_" 
                    + std::to_string(this->collision_id[t_idx][coll_idx]) + ".dat");
                if (!file) {
                    std::cerr << "Error opening file \n";
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                file << "Amount collisions, amount collidable,  Energy loss (J), Incident Energy (J) \n";
                file << std::scientific << std::setprecision(8);
                file << this->total_amount_collisions[t_idx][coll_idx] << "\t" << this->total_amount_collidable_particles << "\t"
                << this->total_energy_loss[t_idx][coll_idx] * 0.5 << "\t"
                << this->total_incident_energy[t_idx][coll_idx] * 0.5 * primary_particle.mass << "\n";
                file.close();
            }     
        }
    }
}


void null_collider::gather_mpi() {
    MPI_Allreduce(MPI_IN_PLACE, &this->total_amount_collidable_particles, 1, mpi_vars::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    for (int t_idx = 0; t_idx < this->number_targets; t_idx++) {
        std::fill(this->total_incident_energy[t_idx].begin(), this->total_incident_energy[t_idx].end(), 0);
        std::fill(this->total_energy_loss[t_idx].begin(), this->total_energy_loss[t_idx].end(), 0);
        std::fill(this->total_amount_collisions[t_idx].begin(), this->total_amount_collisions[t_idx].end(), 0);
        for (int coll_idx = 0; coll_idx < this->number_collisions_per_target[t_idx]; coll_idx++) {
            for (int i_thread = 0; i_thread < omp_get_max_threads(); i_thread++) {
                this->total_incident_energy[t_idx][coll_idx] += this->total_incident_energy_thread[i_thread][t_idx][coll_idx];
                this->total_energy_loss[t_idx][coll_idx] += this->total_energy_loss_thread[i_thread][t_idx][coll_idx];
                this->total_amount_collisions[t_idx][coll_idx] += this->total_amount_collisions_thread[i_thread][t_idx][coll_idx];
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, this->total_incident_energy[t_idx].data(), this->number_collisions_per_target[t_idx], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, this->total_energy_loss[t_idx].data(), this->number_collisions_per_target[t_idx], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, this->total_amount_collisions[t_idx].data(), this->number_collisions_per_target[t_idx], mpi_vars::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    }
    
}

void null_collider::order_collisions() {
    // Reorder collisions based on most collisions to avoid extra looping through collisions
    // do it every diagnostic step
    #pragma omp master 
    {
        if (this->number_targets > 0) {
            std::vector<int> indices_target(this->number_targets);
            std::vector<size_t> total_target_number_collisions(this->number_targets, 0);
            for (int t_idx = 0; t_idx < this->number_targets; t_idx++) {
                std::vector<std::vector<int>> product_indices_copy = this->product_indices[t_idx];
                std::vector<int> collision_type_per_target_copy = this->collision_type_per_target[t_idx];
                std::vector<int> collision_id_copy = this->collision_id[t_idx];
                std::vector<double> energy_threshold_copy = this->energy_threshold[t_idx];
                std::vector<std::vector<double>> sigma_array_copy = this->sigma_array[t_idx];

                indices_target[t_idx] = t_idx;
                std::vector<int> indices_collisions(this->number_collisions_per_target[t_idx]);
                for (int coll_idx=0; coll_idx< this->number_collisions_per_target[t_idx];coll_idx++){
                    indices_collisions[coll_idx] = coll_idx;
                    total_target_number_collisions[t_idx] += this->total_amount_collisions[t_idx][coll_idx];     
                }  
                // sort indicies based on frequency
                std::vector<size_t>& temp_array = this->total_amount_collisions[t_idx];
                std::sort(indices_collisions.begin(), indices_collisions.end(), [&temp_array](int i, int j) {
                    return temp_array[i] > temp_array[j];});

                // sort collisions within target
                for (int coll_idx=0; coll_idx< this->number_collisions_per_target[t_idx];coll_idx++){
                    int idx = indices_collisions[coll_idx];
                    this->collision_id[t_idx][coll_idx] = collision_id_copy[idx];
                    this->collision_type_per_target[t_idx][coll_idx] = collision_type_per_target_copy[idx];
                    this->energy_threshold[t_idx][coll_idx] = energy_threshold_copy[idx];
                    this->product_indices[t_idx][coll_idx] = product_indices_copy[idx];
                    this->sigma_array[t_idx][coll_idx] = sigma_array_copy[idx];
                }   
            }

            std::vector<double> reduced_mass_copy = this->reduced_mass;
            std::vector<double> reduced_mass_ionization_copy = this->reduced_mass_ionization;
            std::vector<int> target_idx_copy = this->target_idx;
            std::vector<int> number_collisions_per_target_copy = this->number_collisions_per_target;
            std::vector<std::vector<std::vector<int>>> product_indices_copy = this->product_indices;
            std::vector<std::vector<int>> collision_type_per_target_copy = this->collision_type_per_target;
            std::vector<std::vector<int>> collision_id_copy = this->collision_id;
            std::vector<std::vector<double>> energy_threshold_copy = this->energy_threshold;
            std::vector<std::vector<std::vector<double>>> sigma_array_copy = this->sigma_array;

            std::sort(indices_target.begin(), indices_target.end(), [&total_target_number_collisions](int i, int j) {
                return total_target_number_collisions[i] > total_target_number_collisions[j];});

            for (int t_idx = 0; t_idx < this->number_targets; t_idx++) {
                int idx = indices_target[t_idx];
                this->number_collisions_per_target[t_idx] = number_collisions_per_target_copy[idx];
                this->target_idx[t_idx] = target_idx_copy[idx];
                this->reduced_mass[t_idx] = reduced_mass_copy[idx];
                this->reduced_mass_ionization[t_idx] = reduced_mass_ionization_copy[idx];
                this->collision_id[t_idx] = collision_id_copy[idx];
                this->collision_type_per_target[t_idx] = collision_type_per_target_copy[idx];
                this->energy_threshold[t_idx] = energy_threshold_copy[idx];
                this->product_indices[t_idx] = product_indices_copy[idx];
                this->sigma_array[t_idx] = sigma_array_copy[idx]; 
            }

        }
    }
}

void null_collider::reset_diagnostics(int thread_id) {
    if (this->number_targets > 0) {
        #pragma omp master 
        {
            this->total_amount_collidable_particles = 0;
        }
        for (int t_idx = 0; t_idx < this->number_targets; t_idx++){
            std::fill(this->total_incident_energy_thread[thread_id][t_idx].begin(), this->total_incident_energy_thread[thread_id][t_idx].end(), 0.0);
            std::fill(this->total_energy_loss_thread[thread_id][t_idx].begin(), this->total_energy_loss_thread[thread_id][t_idx].end(), 0.0);
            std::fill(this->total_amount_collisions_thread[thread_id][t_idx].begin(), this->total_amount_collisions_thread[thread_id][t_idx].end(), 0);
        }
    }
}


std::vector<null_collider> read_null_collision_inputs(const std::string& directory_path, const std::vector<charged_particle> &particle_list, const std::vector<target_particle> &target_particle_list){
    if (mpi_vars::mpi_rank==0){
        std::cout << " "<< std::endl;
        std::cout << "Reading null collision inputs "<< std::endl;
        std::cout << "---------------------------------------- "<< std::endl;
    }

    
    int number_charged_particles = particle_list.size();
    int number_target_particles = target_particle_list.size();
    std::vector<int> number_targets(number_charged_particles, 0);
    std::vector<std::vector<int>> target_indices(number_charged_particles);
    std::vector<std::vector<std::vector<std::vector<int>>>> product_indices(number_charged_particles);
    std::vector<std::vector<int>> number_collisions_per_target(number_charged_particles);
    std::vector<std::vector<std::vector<int>>> collision_type_per_target(number_charged_particles);
    std::vector<std::vector<std::vector<double>>> E_threshold_per_primary_collision(number_charged_particles);
    std::vector<std::vector<std::vector<std::vector<double>>>> energy_arrays(number_charged_particles);
    std::vector<std::vector<std::vector<std::vector<double>>>> sigma_arrays(number_charged_particles);
    double E_threshold, temp_var, E_scaling, sigma_scaling;

    for (int rank_num = 0;rank_num < mpi_vars::mpi_size; rank_num++) {
        if (rank_num == mpi_vars::mpi_rank){
            // Open the directory
            DIR* dir = opendir(directory_path.c_str());
            if (!dir) {
                perror("opendir");
                exit(EXIT_FAILURE);
            }

            struct dirent* entry;
            if ((entry = readdir(dir)) == nullptr && mpi_vars::mpi_rank == 0) {
                std::cout << "No binary collisions included!" << std::endl;
            }
            while ((entry = readdir(dir)) != nullptr) {
                if (entry->d_type == DT_REG) {  // regular file
                    std::string filename = directory_path + entry->d_name;
                    std::ifstream file(filename);
                    if (!file) {
                        std::cerr << "Error: Unable to open file " << filename << std::endl;
                    }
                    std::string line;
                    std::getline(file, line);
                    std::istringstream iss(line);
                    while (line.find("END") == std::string::npos) {
                        if (line.find("REACTION") != std::string::npos) {
                            std::string reaction_string;
                            std::getline(file, line);
                            iss.clear();
                            iss.str(line);
                            iss >> reaction_string;
                            if (mpi_vars::mpi_rank==0) {std::cout << reaction_string << std::endl;}
                            size_t arrow_pos = reaction_string.find("->");
                            std::string reactant_string = reaction_string.substr(0,arrow_pos);
                            std::string product_string = reaction_string.substr(arrow_pos+2);
                            std::regex pattern(R"(\[([^\]]+)\])");
                            std::smatch match;
                            int primary_idx, target_idx; 
                            int i;
                            
                            // Get indices of the primary and target particle
                            while (std::regex_search(reactant_string, match, pattern)){
                                
                                for (i=0;i<number_charged_particles;i++){
                                    if (particle_list[i].name == match[1]){
                                        primary_idx = i;     
                                        break;
                                    }
                                }
                                for (i=0;i<number_target_particles;i++){
                                    if (target_particle_list[i].name == match[1]){
                                        target_idx = i;
                                        break;
                                    }
                                }
                                reactant_string = match.suffix().str();
                            }
                            bool target_exists = false;
                            int target_collision_indx = 0;
                            for (int k = 0; k < number_targets[primary_idx]; k++) {
                                if (target_idx == target_indices[primary_idx][k]) {
                                    target_exists = true;
                                    target_collision_indx = k;
                                    break;
                                }
                            }
                            if (target_exists) {
                                number_collisions_per_target[primary_idx][target_collision_indx]++;
                            } else {
                                number_targets[primary_idx]++;
                                target_indices[primary_idx].push_back(target_idx);
                                number_collisions_per_target[primary_idx].push_back(1);
                            }

                            std::vector<int> local_product_indx;
                            int product_target_idx = -1;
                            while (std::regex_search(product_string, match, pattern)){
                                for (i=0;i<number_charged_particles;i++){
                                    if (particle_list[i].name == match[1]){
                                        local_product_indx.push_back(i);
                                    }
                                }
                                for (i=0;i<number_target_particles;i++){
                                    if (target_particle_list[i].name == match[1]){
                                        product_target_idx = i;
                                    }
                                }
                                product_string = match.suffix().str();
                            }
                            if (product_target_idx > - 1) {
                                local_product_indx.push_back(product_target_idx); // put target index as last
                            }
                            
                            
                            std::getline(file, line);
                            iss.clear();
                            iss.str(line);
                            iss >> E_threshold >> temp_var;
                            std::getline(file, line);
                            iss.clear();
                            iss.str(line);
                            iss >> E_scaling >> sigma_scaling;
                            std::string coll_string;
                            std::getline(file, line);
                            iss.clear();
                            iss.str(line);
                            iss >> coll_string;
                            int collision_type = 0;
                            if (coll_string == "ELASTIC") {
                                collision_type = 1;
                            } else if (coll_string == "EXCITATION"){
                                collision_type = 3;
                            } else if (coll_string == "IONIZATION"){
                                collision_type = 2;
                                // make sure primary particle is first, then followed by by-product
                                for (i=0;i<3;i++){
                                    if (primary_idx == local_product_indx[i]){
                                        // switch primary to front
                                        int temp_idx = local_product_indx[0];
                                        local_product_indx[0] = primary_idx;
                                        local_product_indx[i] = temp_idx;
                                        break;
                                    }
                                }
                                for (i=1;i<3;i++){
                                    if (particle_list[local_product_indx[i]].mass == constants::electron_mass){
                                        // switch electron to back
                                        int local_idx = local_product_indx[i];
                                        int end_idx = local_product_indx[2];
                                        local_product_indx[2] = local_idx;
                                        local_product_indx[i] = end_idx;
                                        break;
                                    }
                                }
                                // Make sure masses of product and electron equal to target mass
                                double target_mass = target_particle_list[target_idx].mass;
                                double product_mass = particle_list[local_product_indx[1]].mass;
                                if (std::abs(product_mass + constants::electron_mass - target_mass)/target_mass > 1e-6) {
                                    std::cout << "Ionization masses do not add up!" << std::endl;
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                local_product_indx.pop_back();
                            } else if (coll_string == "CHARGEEXCHANGE"){
                                collision_type = 4;
                            }
                            if (target_exists) {
                                collision_type_per_target[primary_idx][target_collision_indx].push_back(collision_type);
                                E_threshold_per_primary_collision[primary_idx][target_collision_indx].push_back(E_threshold);
                                product_indices[primary_idx][target_collision_indx].push_back(local_product_indx);
                            } else {
                                collision_type_per_target[primary_idx].push_back({collision_type});
                                E_threshold_per_primary_collision[primary_idx].push_back({E_threshold});
                                product_indices[primary_idx].push_back({local_product_indx});
                            }
                            
                            
                            std::getline(file, line);
                            while (line.find("---------") == std::string::npos) {
                                std::getline(file, line);
                            }
                            double energy_point, sigma_point;
                            std::vector<double> local_energy_array;
                            std::vector<double> local_sigma_array;
                            std::getline(file, line);
                            while (line.find("---------") == std::string::npos) {
                                iss.clear();
                                iss.str(line);
                                iss >> energy_point >> sigma_point;
                                local_energy_array.push_back(energy_point * E_scaling);
                                local_sigma_array.push_back(sigma_point * sigma_scaling);
                                std::getline(file, line);
                            }
                            // find which particle indices they match to
                            if (target_exists) {
                                energy_arrays[primary_idx][target_collision_indx].push_back(local_energy_array);
                                sigma_arrays[primary_idx][target_collision_indx].push_back(local_sigma_array);
                            } else {
                                energy_arrays[primary_idx].push_back({local_energy_array});
                                sigma_arrays[primary_idx].push_back({local_sigma_array});
                            }
                
                        }

                        // go to next line
                        std::getline(file, line);
                    }
                    file.close();
                    
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    

    // Concatenate into final energy/sigma arrays
    // Overall seems likely better 1 binary search per particle regardless of how big array becomes because of O(log(n)) time
    std::vector<int> number_energy_points(number_charged_particles);
    std::vector<std::vector<double>> total_energy_array(number_charged_particles);
    std::vector<std::vector<std::vector<std::vector<double>>>> total_sigma_array(number_charged_particles);
    for (int k=0; k<number_charged_particles; k++){
        
        if (number_targets[k] > 0) {
            total_sigma_array[k].resize(number_targets[k]);
            std::vector<std::vector<int>> current_indx_collision_array(number_targets[k]);
            for (int t_idx = 0; t_idx < number_targets[k]; t_idx++){
                current_indx_collision_array[t_idx].resize(number_collisions_per_target[k][t_idx], 0);
                total_sigma_array[k][t_idx].resize(number_collisions_per_target[k][t_idx]);
            }
            
            bool not_finished = true;
            double curr_min_value;
            // Find current minimum energy value
            // concatenate energy array
            int t_idx_temp, coll_idx_temp;
            while (not_finished) {
                not_finished = false;
                for (int t_idx = 0; t_idx < number_targets[k]; t_idx++) {
                    for (int coll_idx=0; coll_idx < number_collisions_per_target[k][t_idx]; coll_idx++){
                        int array_size = energy_arrays[k][t_idx][coll_idx].size();
                        int current_idx = current_indx_collision_array[t_idx][coll_idx];
                        double current_energy = energy_arrays[k][t_idx][coll_idx][current_idx];
                        if (current_idx < array_size) {
                            if (!not_finished) {
                                // first usable index is min value
                                curr_min_value = current_energy;
                                t_idx_temp = t_idx;
                                coll_idx_temp = coll_idx;
                                not_finished = true;
                            } else {
                                if (current_energy < curr_min_value){
                                    curr_min_value = current_energy;
                                    t_idx_temp = t_idx;
                                    coll_idx_temp = coll_idx;
                                } else if (current_energy == curr_min_value) {
                                    // if equal, since we already have a value we can use, we just increase the index for that collision array
                                    current_indx_collision_array[t_idx][coll_idx]++;
                                }
                            }
                        } 
                    }
                } 
                if (not_finished) {
                    total_energy_array[k].push_back(curr_min_value);
                    current_indx_collision_array[t_idx_temp][coll_idx_temp]++;
                }
            }
            
            // Now interpolate to sigma array for each collision
            for (int t_idx = 0; t_idx < number_targets[k]; t_idx++) {
                for (int coll_idx=0; coll_idx < number_collisions_per_target[k][t_idx]; coll_idx++){
                    int lower_idx = 0;
                    double interp_d, temp_var;
                    total_sigma_array[k][t_idx][coll_idx].resize(total_energy_array[k].size());
                    for (int u = 0; u < total_energy_array[k].size();u++){
                        curr_min_value = total_energy_array[k][u];
                        if (curr_min_value < E_threshold_per_primary_collision[k][t_idx][coll_idx]) {
                            // outside of lower energy array
                            total_sigma_array[k][t_idx][coll_idx][u] = 0.0;
                        } else if (curr_min_value < energy_arrays[k][t_idx][coll_idx].back()) {
                            while (energy_arrays[k][t_idx][coll_idx][lower_idx] <= curr_min_value){
                                lower_idx++;
                            }
                            lower_idx--;
                            if (energy_arrays[k][t_idx][coll_idx][lower_idx] == curr_min_value){
                                total_sigma_array[k][t_idx][coll_idx][u] = sigma_arrays[k][t_idx][coll_idx][lower_idx];
                            } else{
                                temp_var = curr_min_value - energy_arrays[k][t_idx][coll_idx][lower_idx];
                                interp_d = temp_var / (energy_arrays[k][t_idx][coll_idx][lower_idx+1] - energy_arrays[k][t_idx][coll_idx][lower_idx]);
                                // linear interpolate sigma
                                total_sigma_array[k][t_idx][coll_idx][u] = sigma_arrays[k][t_idx][coll_idx][lower_idx] * (1.0 - interp_d) + sigma_arrays[k][t_idx][coll_idx][lower_idx+1] * (interp_d);
                            }
                        } else {
                            // outside of max energy, set to maximum sigma
                            total_sigma_array[k][t_idx][coll_idx][u] = sigma_arrays[k][t_idx][coll_idx].back();
                        }
                    }
                }
            }
            number_energy_points[k] = total_energy_array[k].size();
        }
    }

    // Sort collisions on most likely collisions based on sigma * v
    // Just create new in memory, will be deallocated
    
    for (int k=0; k< number_charged_particles; k++){
        if (number_targets[k] > 0) {
            // copy vectors for sorting after
            std::vector<double> freq_max_target(number_targets[k], 0.0);
            std::vector<int> indices_target(number_targets[k]);
            for (int t_idx = 0; t_idx < number_targets[k]; t_idx++){
                std::vector<std::vector<int>> product_indices_copy = product_indices[k][t_idx];
                std::vector<int> collision_type_per_target_copy = collision_type_per_target[k][t_idx];
                std::vector<double> E_threshold_per_primary_collision_copy = E_threshold_per_primary_collision[k][t_idx];
                std::vector<std::vector<double>> total_sigma_array_copy = total_sigma_array[k][t_idx];
                indices_target[t_idx] = t_idx;
                freq_max_target[t_idx] = 0.0;
                std::vector<double> freq_max_coll(number_collisions_per_target[k][t_idx], 0);
                std::vector<int> indices_coll(number_collisions_per_target[k][t_idx]);
                double density = target_particle_list[target_indices[k][t_idx]].average_density;
                double target_mass = target_particle_list[target_indices[k][t_idx]].mass;
                if (particle_list[k].number_velocity_coordinates < 3) {
                    std::cout << "Primary particle " << particle_list[k].name << " does not have 3D velocity!" << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                double primary_mass = particle_list[k].mass;
                double red_mass = (primary_mass*target_mass)/(primary_mass + target_mass);
                for (int coll_idx = 0; coll_idx < number_collisions_per_target[k][t_idx]; coll_idx++){
                    double v_r;
                    indices_coll[coll_idx] = coll_idx;
                    freq_max_coll[coll_idx] = 0.0;
                    for (int p=0;p<total_energy_array[k].size();p++) {
                        double energy_cm = total_energy_array[k][p];
                        double sigma = total_sigma_array[k][t_idx][coll_idx][p];
                        v_r = std::sqrt(energy_cm * 2.0 * constants::elementary_charge / red_mass); // convert to relative velocity
                        freq_max_coll[coll_idx] = std::max(freq_max_coll[coll_idx], density * v_r * sigma); // get max for this collision
                    }
                    freq_max_target[t_idx] += freq_max_coll[coll_idx];
                }
                // sort indicies based on frequency
                std::sort(indices_coll.begin(), indices_coll.end(), [&freq_max_coll](int i, int j) {
                    return freq_max_coll[i] > freq_max_coll[j];});

                // sort collisions within target
                for (int coll_idx = 0; coll_idx < number_collisions_per_target[k][t_idx]; coll_idx++){
                    int idx = indices_coll[coll_idx];
                    product_indices[k][t_idx][coll_idx] = product_indices_copy[idx];
                    collision_type_per_target[k][t_idx][coll_idx] = collision_type_per_target_copy[idx];
                    E_threshold_per_primary_collision[k][t_idx][coll_idx] = E_threshold_per_primary_collision_copy[idx];
                    total_sigma_array[k][t_idx][coll_idx] = total_sigma_array_copy[idx];
                }
                
            }
            
            // Create copies to then switch target order by maximum amount of collisions
            std::vector<int> target_indices_copy = target_indices[k];
            std::vector<int> number_collisions_per_target_copy = number_collisions_per_target[k];
            std::vector<std::vector<std::vector<int>>> product_indices_copy = product_indices[k];
            std::vector<std::vector<int>> collision_type_per_target_copy = collision_type_per_target[k];
            std::vector<std::vector<double>> E_threshold_per_primary_collision_copy = E_threshold_per_primary_collision[k];
            std::vector<std::vector<std::vector<double>>> total_sigma_array_copy = total_sigma_array[k];

            // sort indicies based on frequency
            std::sort(indices_target.begin(), indices_target.end(), [&freq_max_target](int i, int j) {
                return freq_max_target[i] > freq_max_target[j];});

            // sort collisions within particle
            for (int t_idx = 0; t_idx < number_targets[k]; t_idx++){
                int idx = indices_target[t_idx];
                target_indices[k][t_idx] = target_indices_copy[idx];
                number_collisions_per_target[k][t_idx] = number_collisions_per_target_copy[idx];
                product_indices[k][t_idx] = product_indices_copy[idx];
                collision_type_per_target[k][t_idx] = collision_type_per_target_copy[idx];
                E_threshold_per_primary_collision[k][t_idx] = E_threshold_per_primary_collision_copy[idx];
                total_sigma_array[k][t_idx] = total_sigma_array_copy[idx];
            }
            
        }
    }


    std::vector<null_collider> binary_collision_list;
    binary_collision_list.reserve(number_charged_particles);
    for (int k=0; k<number_charged_particles; k++){
        std::vector<double> red_mass(number_targets[k]), red_mass_ion(number_targets[k]);
        double primary_mass = particle_list[k].mass;
        for (int i = 0; i < number_targets[k];i++){
            double target_mass = target_particle_list[target_indices[k][i]].mass;
            red_mass[i] = (target_mass * primary_mass)/(target_mass + primary_mass);
            red_mass_ion[i] = 1.0 / (1.0/primary_mass + 1.0/(target_mass - constants::electron_mass) + 1.0/constants::electron_mass);
        }
        binary_collision_list.emplace_back(k, number_targets[k], target_indices[k], number_collisions_per_target[k], 
            total_sigma_array[k], total_energy_array[k], E_threshold_per_primary_collision[k],
            collision_type_per_target[k], product_indices[k], red_mass, red_mass_ion);
    }
    
    return binary_collision_list;

}