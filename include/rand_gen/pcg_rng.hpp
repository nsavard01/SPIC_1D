#pragma once
#include <cstdint>
// Define multiplier for conversion to [0,1] double


// Declare PRNG state as a global variable
// extern uint64_t pcg_state;

// Function to initialize the PRNG state for each thread.  seed_stream = 0 draws a fresh
// non-reproducible seed; seed_stream >= 1 selects a reproducible, independent stream, so
// that repeats of the same case differ only in the random sequence.
void initialize_pcg(int seed_stream);

// Function to generate a random number in [0,1]
double pcg32_random_r();

uint64_t get_pcg_state();


