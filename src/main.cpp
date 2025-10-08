#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <fstream>
#include <sstream>
#include "globals/mpi_vars.hpp"
#include "simulation/simulation.hpp"
#include <mkl.h>


int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_vars::mpi_size);  // Get total processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_vars::mpi_rank);  // Get current rank

    // Determine the MPI datatype for size_t based on its size
    if (sizeof(size_t) == sizeof(unsigned int)) {
        mpi_vars::mpi_size_t_type = MPI_UNSIGNED;
    } else if (sizeof(size_t) == sizeof(unsigned long)) {
        mpi_vars::mpi_size_t_type = MPI_UNSIGNED_LONG;
    } else if (sizeof(size_t) == sizeof(unsigned long long)) {
        mpi_vars::mpi_size_t_type = MPI_UNSIGNED_LONG_LONG;
    } else {
        std::cout << "Unsupported size_t type!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    
    simulation simulator;
    simulator.setup();
    simulator.run();
    simulator.averaging();
    MPI_Finalize();
    return 0;
}