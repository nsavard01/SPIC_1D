
#include "particles/target_particle.hpp"
#include "globals/constants.hpp"
#include "globals/mpi_vars.hpp"
#include <cmath>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <mpi.h>
#include <numeric>
#include <algorithm>

target_particle::target_particle(double mass_in, double temp_in, double density_in, double v_drift_in, std::string name_in)
    {
    this->name = name_in;
    this->mass = mass_in;
    this->average_density = density_in;
    this->average_temperature = temp_in;
    this->v_drift = v_drift_in;
    this->v_therm = std::sqrt(this->average_temperature * constants::k_boltz / this->mass);
    
}


void target_particle::print_out() const {
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "Target particle name: " << this->name << std::endl;
        std::cout << "Mass (kg): " << this->mass << std::endl;
        std::cout << "Average density (m^-3): " << this->average_density << std::endl;
        std::cout << "Average temperature (K): " << this->average_temperature << std::endl;
        std::cout << " " << std::endl;
    }
}

void target_particle::initialize_diagnostic_files(const std::string& dir_name) const {
    if (mpi_vars::mpi_rank == 0) {
        std::ofstream file(dir_name + "/target_particles/" + this->name + "/particle_properties.dat");

        // Write header (optional)
        file << "Particle Symbol, Particle Mass (kg) \n";
        file << std::scientific << std::setprecision(8);
        file << this->name << "\t"
            << this->mass << "\t"
            <<"\n";


        file.close();

        file.open(dir_name + "/target_particles/" + this->name + "/ave_diagnostics.dat");
        if (!file) {
            std::cerr << "Error opening file for energy particle \n";
            return;
        }

        file << "density (1/m^3), temperature (K) \n";

        file.close();
    }
}

void target_particle::write_diagnostics(const std::string& dir_name, int diag_number) const {
    if (mpi_vars::mpi_rank == 0) {
    
        std::ofstream file(dir_name + "/target_particles/" + this->name + "/ave_diagnostics.dat", std::ios::app);
        if (!file) {
            std::cerr << "Error opening file for target particle \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        file << std::scientific << std::setprecision(8);
        file << this->average_density << "\t"
        << this->average_temperature << "\n";

        file.close();
    }
}




std::vector<target_particle> read_target_particle_inputs(const std::string& directory_path, const domain& world){
    std::vector<target_particle> target_particle_list;
    std::vector<double> mass_in, n_ave, temp_in, v_drift;
    std::vector<std::string> particle_names;
    std::vector<int> index_order;
    int count_number_particles = 0;
    if (mpi_vars::mpi_rank == 0) {
        std::cout << " "<< std::endl;
        std::cout << "Reading target particle inputs "<< std::endl;
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
                    double mass;
                    iss >> mass;
                    mass = mass * constants::mass_amu;
                    mass_in.push_back(mass);
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
                    double drift;
                    iss >> drift;
                    v_drift.push_back(drift);
                    file.close();
                }
            }   
        
            closedir(dir);

            // Initialize index_order with indices [0, 1, 2, ..., count_number_particles - 1]
            index_order.resize(count_number_particles);
            std::iota(index_order.begin(), index_order.end(), 0);

            
            // Sort indices based on initial density
            std::sort(index_order.begin(), index_order.end(), [&](int i, int j) {
                return (n_ave[i]) > (n_ave[j]);
            });
        
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    

    //Create the particles in the sorted order 
    for (int num_part = 0; num_part < count_number_particles; num_part++){
        int i = index_order[num_part];
        target_particle temp_particle(mass_in[i], temp_in[i], n_ave[i], v_drift[i], particle_names[i]);
        target_particle_list.push_back(temp_particle);
    }
    
    for (int i = 0; i < target_particle_list.size(); i++) {
        target_particle_list[i].print_out();
    }

    return target_particle_list;

}


