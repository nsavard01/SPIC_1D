#pragma once
#include <cstdint>
// Define multiplier for conversion to [0,1] double


// Declare PRNG state as a global variable
// extern uint64_t pcg_state;

// Function to initialize the PRNG state for each thread
void initialize_pcg(bool pre_determined);

// Function to generate a random number in [0,1]
double pcg32_random_r();

uint64_t get_pcg_state();


