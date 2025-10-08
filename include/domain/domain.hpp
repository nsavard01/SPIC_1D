// include/domain.hpp
#pragma once
#include <vector>
#include <string>
#include <memory>

class domain {
    
public:

    int number_cells; // number of cells in the domain
    int domain_type;
    int number_nodes; // number of nodes in the domain
    double length_domain; // length of the domain
    double min_dx; // minimum cell size
    std::vector<double> dx_dxi; // cell sizes in the domain
    std::vector<double> grid_nodes; // grid nodes
    std::vector<double> cell_centers; // cell centers
    int right_boundary_condition; // right boundary condition type
    int left_boundary_condition; // left boundary condition type
    virtual ~domain() = default;

    


    virtual void print_out() = 0;
    virtual void write_domain(const std::string& filename) = 0;
    // virtual void read_domain_from_file(const std::string& filename);
};

std::unique_ptr<domain> create_domain_from_file(const std::string& filename, int scheme_type);
