
#pragma once
#include "charged_particle_operator.hpp"
#include "rand_gen/pcg_rng.hpp"
class charged_particle_wall_injector : public charged_particle_operator {

public:
    
    
    std::vector<double> current_density, particle_location, v_therm;
    std::vector<std::vector<double>> v_3D;
    std::vector<std::vector<std::vector<size_t>>> accumulated_number_particles;
    charged_particle_wall_injector(std::vector<charged_particle>& particle_list, std::vector<int>& particle_indx, std::vector<double>& current_density, std::vector<std::vector<double>>& v_3D, 
        std::vector<double>& v_therm, std::vector<double>& particle_location);
    void setup_diagnostics(const std::string& dir_name, std::vector<charged_particle>& particle_list) override;
    void write_diagnostics(const std::string& dir_name, const std::vector<charged_particle>& particle_list) override;
    void reset_diagnostics() override;
    void write_average_diagnostics(const std::string& dir_name, const std::vector<charged_particle>& particle_list) override;
    void print_out() override;
    void run(const int thread_id, const double current_time, const double del_t, std::vector<charged_particle>& particle_list, const domain& world) override;

};


