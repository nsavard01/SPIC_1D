// include/domain.hpp
#pragma once
#include "domain/domain.hpp"
#include <vector>
#include <iostream>

class non_uniform_domain : public domain {
protected:
    bool curvilinear; // curvilinear flag

public:

    
    non_uniform_domain(int num_cells, double length_domain, int left_boundary_condition, int right_boundary_condition,
        int type, int temp_int_1, int temp_int_2, double temp_double_1, double temp_double_2);
    void print_out() override;
    void write_domain(const std::string& filename) override;
    

};

// 
