// include/domain.hpp
#pragma once
#include "domain/domain.hpp"
#include <vector>
#include <iostream>

class uniform_domain : public domain {

public:

    
    uniform_domain(int num_cells, double length_domain, int left_boundary_condition, int right_boundary_condition);

    void print_out() override;
    void write_domain(const std::string& filename) override;


};
