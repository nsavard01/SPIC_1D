#pragma once
#include "rand_gen/pcg_rng.hpp"
#include <cmath>
// function to generate Maxwellian distributions
inline void maxwellian_1D(double& v_x, const double& v_therm, const double& v_drift) {
    double R_1 = pcg32_random_r();
    double R_2 = pcg32_random_r();
    v_x = v_drift + v_therm * std::sqrt(-2.0 * std::log(R_1)) * std::cos(2.0 * M_PI * R_2); 
}

inline void maxwellian_2D(double& v_x, double& v_y, const double& v_therm, const double& v_drift) {
    double R_1 = pcg32_random_r();
    double R_2 = pcg32_random_r();
    double coeff = std::sqrt(-2.0 * std::log(R_1));
    v_x = v_drift + v_therm * coeff * std::cos(2.0 * M_PI * R_2);   
    v_y = v_therm * coeff * std::sin(2.0 * M_PI * R_2);
}

inline void maxwellian_3D(double& v_x, double& v_y, double& v_z, const double& v_therm, const double& v_drift) {
    double R_1 = pcg32_random_r();
    double R_2 = pcg32_random_r();
    double R_3 = pcg32_random_r();
    double R_4 = pcg32_random_r();
    double coeff = std::sqrt(-2.0 * std::log(R_1));
    v_x = v_drift + v_therm * coeff * std::cos(2.0 * M_PI * R_2); 
    v_y = v_therm * coeff * std::sin(2.0 * M_PI * R_2);
    v_z = v_therm * std::sqrt(-2.0 * std::log(R_3)) * std::cos(2.0 * M_PI * R_4);
}

inline void maxwellian_3D_flux(double& v_x, double& v_y, double& v_z, const double& v_therm, const double& v_drift) {
    double R_1 = pcg32_random_r();
    double R_2 = pcg32_random_r();
    double R_3 = pcg32_random_r();
    double R_4 = pcg32_random_r();
    double coeff = std::sqrt(-2.0 * std::log(R_1));
    v_x = v_drift + v_therm * coeff; 
    v_y = v_therm * coeff * std::sin(2.0 * M_PI * R_2);
    v_z = v_therm * std::sqrt(-2.0 * std::log(R_3)) * std::cos(2.0 * M_PI * R_4);
}



