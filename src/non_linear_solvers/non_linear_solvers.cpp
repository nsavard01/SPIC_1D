#include "non_linear_solvers/non_linear_solver.hpp"
#include "globals/mpi_vars.hpp"
#include <iostream>
#include <fstream>
#include <sstream> 
#include <iomanip>

void non_linear_solver::initialize_diagnostic_files(const std::string& filename) const {
    if (mpi_vars::mpi_rank == 0) {
        std::ofstream file(filename + "non_linear_solver_properties.dat");
        if (!file) {
            std::cerr << "Error opening file for non_linear_properties \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Write header (optional)
        file << "absolute tolerance, relative tolerance, max iterations, number unknowns \n";

        file << std::scientific << std::setprecision(8);
        file << this->eps_a << "\t"
        << this->eps_r << "\t"
        << this->max_iterations << "\t"
        << this->number_unknowns <<
        "\n";

        file.close();

        file.open(filename + "/non_linear_solver_diagnostics.dat");
        if (!file) {
            std::cerr << "Error opening file for domain \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Write header (optional)
        file << "solver time (s), residual norm, number iterations \n";

        file.close();
    }
}


void non_linear_solver::write_diagnostics(const std::string& filename) const {
    if (mpi_vars::mpi_rank == 0) {
        std::ofstream file(filename + "/non_linear_solver_diagnostics.dat", std::ios::app);
        

        if (!file) {
            std::cerr << "Error opening file for domain \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Write header (optional)
        file << std::scientific << std::setprecision(8);
        file << this->solver_time << "\t"
        << this->accum_residual_norm << "\t"
        << this->accum_iterations_count <<
        "\n";

        file.close();
    }
}

void read_non_linear_solver_inputs(const std::string& filename, std::vector<int>& int_params, std::vector<double>& double_params){
    if (mpi_vars::mpi_rank == 0) {
        std::cout << " "  << std::endl;
        std::cout << "Reading non-linear inputs: "  << std::endl;
        std::cout << "-------------------------- "  << std::endl;
        std::string line;
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        
        std::getline(file, line);
        std::istringstream iss(line);
        iss >> int_params[0];
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> int_params[1];
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> int_params[2];
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> double_params[0];
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> double_params[1];
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> double_params[2];
        iss.clear();
        file.close();
        std::cout << " "  << std::endl;
    }
}


