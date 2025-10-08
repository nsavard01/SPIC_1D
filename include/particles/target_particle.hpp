#pragma once

#include <vector>
#include <cstddef>
#include <string>
#include "domain/domain.hpp"

class target_particle {
public:
    double mass, average_density, average_temperature, v_therm, accum_energy_change, diffusion_coeff, v_drift;
    std::string name;
    // Accessor for 3D indexing
    target_particle(double mass_in, double temp_in, double density_in, double v_drift_in, std::string name_in);
    void print_out() const;
    void initialize_diagnostic_files(const std::string& dir_name) const;
    void write_diagnostics(const std::string& dir_name, int diag_number) const;
};
std::vector<target_particle> read_target_particle_inputs(const std::string& filename, const domain& world);

