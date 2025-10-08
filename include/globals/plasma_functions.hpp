#pragma once
#include "globals/constants.hpp"
#include <cmath>

// functions for plasma physics simulations

double get_debye_length(double temperature, double density, double charge, double mass) {
    // Calculate the Debye length based on the given parameters
    // temperature in eV, density in m^-3, charge in C, mass in kg
    return std::sqrt(constants::epsilon_0 * constants::elementary_charge * temperature / (density * charge * charge));
}

double get_plasma_frequency(double density) {
    // Calculate the plasma frequency based on the given parameters
    // temperature in eV, density in m^-3
    return std::sqrt(density * constants::elementary_charge * constants::elementary_charge / (constants::epsilon_0 * constants::electron_mass));
}