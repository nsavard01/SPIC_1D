#include "non_linear_solvers/AA_solver.hpp"
#include "globals/mpi_vars.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream> 



AA_solver::AA_solver(double beta, double eps_r, double eps_a, int max_iterations, int m_anderson, int number_unknowns){
    this->beta = beta;
    this->eps_r = eps_r;
    this->eps_a = eps_a;
    this->max_iterations = max_iterations;
    this->m_anderson = m_anderson;
    this->number_unknowns = number_unknowns;
    this->norm_residual.resize(this->m_anderson+1, 0.0);
    this->min_matrix.resize(this->number_unknowns * this->m_anderson);
    this->residual_k.resize(this->m_anderson + 1);
    this->x_k.resize(this->m_anderson + 1);
    for (int i = 0; i < this->m_anderson + 1; i++) {
        this->residual_k[i].resize(this->number_unknowns, 0.0);
        this->x_k[i].resize(this->number_unknowns, 0.0);
    }
    this->accum_iterations_count = 0;
    this->accum_residual_norm = 0.0;
    this->solver_time = 0.0;
}

void AA_solver::initialize_diagnostic_files(const std::string& filename) const {
    if (mpi_vars::mpi_rank == 0) {
        std::ofstream file(filename + "/non_linear_solver_properties.dat");
        if (!file) {
            std::cerr << "Error opening file for non_linear_properties \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Write header (optional)
        file << "solver type, absolute tolerance, relative tolerance, max iterations, number unknowns, m_anderson, beta \n";

        file << std::scientific << std::setprecision(8);
        file << "AA" << "\t" << this->eps_a << "\t"
        << this->eps_r << "\t"
        << this->max_iterations << "\t"
        << this->number_unknowns << "\t"
        << this->m_anderson << "\t"
        << this->beta <<
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

void AA_solver::print_out() const {
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "AA_solver: " << std::endl;
        std::cout << "-------------------------- " << std::endl;
        std::cout << "Beta: " << this->beta << std::endl;
        std::cout << "Epsilon_r: " << this->eps_r << std::endl;
        std::cout << "Epsilon_a: " << this->eps_a << std::endl;
        std::cout << "Max iterations: " << this->max_iterations << std::endl;
        std::cout << "M Anderson: " << this->m_anderson << std::endl;
        std::cout << "Number of unknowns: " << this->number_unknowns << std::endl;
        std::cout << "-------------------------- " << std::endl;
    }
}

