#pragma once
#include <mpi.h>
#include <iostream>

// global physical constants
namespace mpi_vars {
    extern int mpi_rank, mpi_size;
    extern MPI_Datatype mpi_size_t_type;
}