#include "rand_gen/pcg_rng.hpp"
#include "globals/mpi_vars.hpp"
#include <omp.h>
#include <random>
#include <stdio.h>
#include <iostream>

// Define PRNG state globally, but make it thread-private
static const double mult_factor_PCG = 1.0 / (static_cast<double>(UINT32_MAX) + 2.0);
static const uint64_t multiplier = 6364136223846793005ULL;
static const uint64_t increment = 1442695040888963407ULL;
static uint64_t pcg_state;
#pragma omp threadprivate(pcg_state)

void initialize_pcg(bool pre_determined) {
    uint64_t seed;
    int thread_id;
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "----------------------------------------" << std::endl;
        if (pre_determined) {
            std::cout << "Using pre-determined seed for PCG RNG." << std::endl;
        } else {
            std::cout << "Using random seed for PCG RNG." << std::endl;
        }
        std::cout << "----------------------------------------" << std::endl;
    }
    #pragma omp parallel private(seed, thread_id) 
    {
    
        thread_id = omp_get_thread_num();
        if (pre_determined) {
            seed = thread_id*multiplier + increment;
        }
        else {
            thread_local std::random_device rd;
            thread_local std::mt19937_64 rng(rd());
            std::uniform_int_distribution<uint64_t> dist(0, UINT64_MAX);
            seed = dist(rng);
        }
        
        pcg_state = seed; // Unique seed per thread
    }
}

double pcg32_random_r()
{   
    // *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
    // Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
    // Convert state to uint64_t
    uint64_t oldstate = pcg_state;
    // Advance internal state
    pcg_state = oldstate * multiplier + increment;
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = static_cast<uint32_t>(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = static_cast<uint32_t>(oldstate >> 59u);
    uint32_t res = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    return ((static_cast<double>(res)) + 1.0) * mult_factor_PCG;
}

uint64_t get_pcg_state() {
    return pcg_state;
}
