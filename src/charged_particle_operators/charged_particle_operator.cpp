
#include "charged_particle_operators/charged_particle_operator.hpp"
#include "charged_particle_operators/charged_particle_wall_injector.hpp"
#include <iomanip>
#include <dirent.h>
#include <fstream>
#include <sstream>


std::vector<std::unique_ptr<charged_particle_operator>> read_particle_operators(const std::string& directory_path, std::vector<charged_particle>& particle_list, const domain& world) {
    std::vector<std::unique_ptr<charged_particle_operator>> output;
    if (mpi_vars::mpi_rank == 0) {
        std::cout << " "<< std::endl;
        std::cout << "Reading charged particle operations "<< std::endl;
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
                    std::string filename = directory_path + entry->d_name;
                    if (!filename.compare(directory_path + "wall_injection.inp")) {
                        std::vector<double> current_density, v_therm, particle_location;
                        std::vector<int> particle_indx;
                        std::vector<std::vector<double>> v_3D;
                        std::string filename = directory_path + entry->d_name;
                        std::string line;
                        std::ifstream file(filename);
                        if (!file) {
                            std::cerr << "Error: Unable to open file " << filename << std::endl;
                            exit(EXIT_FAILURE);
                        }
                        std::getline(file, line);
                        std::istringstream iss(line);
                        while (line.find("END") == std::string::npos) {
                            if (line.find("----") != std::string::npos) {
                                std::vector<double> v_local;
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                std::string name;
                                iss >> name;
                                int indx = -1;
                                for (int i = 0; i < particle_list.size(); i++) {
                                    if (particle_list[i].name == name) {
                                        indx = i;
                                        break;
                                    }
                                }
                                if (indx < 0) {
                                    std::cout << "particle " << name << " does not exist with index " << indx << std::endl;
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                particle_indx.push_back(indx);
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                int node;
                                iss >> node;
                                if (node != 0 && node != world.number_nodes) {
                                    std::cout << "ERROR: Node for particle injection not on boundary!" << std::endl;
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                } else if (node == 0) {
                                    if ((world.left_boundary_condition != 1) && (world.left_boundary_condition != 4)){
                                        std::cout << "WARNING: Leftmost node for particle injection not on metallic boundary!" << std::endl;
                                    }
                                } else if (node == world.number_nodes) {
                                    if ((world.right_boundary_condition != 1) && (world.right_boundary_condition != 4)){
                                        std::cout << "WARNING: Rightmost node for particle injection not on metallic boundary!" << std::endl;
                                    }
                                }
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                double J;
                                iss >> J;
                                J = std::abs(J);
                                current_density.push_back(J);
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                double v_x;
                                iss >> v_x;
                                if (v_x > 0 && node == world.number_nodes) {
                                    std::cout << "ERROR: V_x > 0 put at rightmost node for particle injection!" << std::endl;
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                } else if (v_x < 0 && node == 0) {
                                    std::cout << "ERROR: V_x < 0 put at leftmost node for particle injection!" << std::endl;
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                // if v_x == 0, then make tiny number so direction is known
                                if (v_x == 0.0 && node == world.number_nodes){
                                    v_x = - constants::machine_eps;
                                } else if (v_x == 0.0 && node == 0){
                                    v_x = constants::machine_eps;
                                }
                                v_local.push_back(v_x);
                                // To avoid issues with other particle operations near boundary, put slightly within domain
                                double new_position;
                                if (v_x > 0) {
                                    new_position = std::nextafter(double(node), node+1);
                                } else {
                                    new_position = std::nextafter(double(node), node-1);
                                }
                                particle_location.push_back(new_position);
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                double v_y;
                                iss >> v_y;
                                v_local.push_back(v_y);
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                double v_z;
                                iss >> v_z;
                                v_local.push_back(v_z);
                                v_3D.push_back(v_local);
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                double temperature;
                                iss >> temperature;
                                double v_therm_local = std::sqrt(temperature * std::abs(particle_list[indx].q_over_m));
                                v_therm.push_back(v_therm_local);
                                iss.clear();
                                std::getline(file, line);
                            }
                            
                            std::getline(file, line);
                        }
                        file.close();
                        output.push_back(std::make_unique<charged_particle_wall_injector>(particle_list, particle_indx, current_density, v_3D, v_therm, particle_location));
                    }
                }
            }   
        
            closedir(dir);
        
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    return output;

}
