#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include "globals/write_functions.hpp"
#include "globals/mpi_vars.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <memory>
#include <iomanip>

// Uniform domain constructor

uniform_domain::uniform_domain(int num_cells, double length_domain, int left_boundary_condition, int right_boundary_condition) {
    this->domain_type = 0;
    this->number_cells = num_cells;
    this->length_domain = length_domain;
    this->left_boundary_condition = left_boundary_condition;
    this->right_boundary_condition = right_boundary_condition;
    this->min_dx = length_domain / num_cells;
    this->number_nodes = num_cells + 1;
    this->grid_nodes.resize(this->number_nodes);
    this->cell_centers.resize(this->number_cells);
    for (size_t i = 0; i < this->number_nodes; i++) {
        this->grid_nodes[i] = i * this->min_dx;
    }
    for (size_t i = 0; i < this->number_cells; i++) {
        this->cell_centers[i] = (this->grid_nodes[i] + this->grid_nodes[i + 1]) * 0.5;
    }
}

void uniform_domain::print_out() {
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "Uniform domain with " << this->number_cells << " cells and " << this->number_nodes << " nodes." << std::endl;
        std::cout << "Length of the domain: " << this->length_domain << std::endl;
        std::cout << "Minimum cell size: " << this->min_dx << std::endl;
        std::cout << "Left boundary condition: " << this->left_boundary_condition << std::endl;
        std::cout << "Right boundary condition: " << this->right_boundary_condition << std::endl;
        std::cout << "--------------------------------" << std::endl;
    }
}

// Non-uniform domain constructor

double sinusoidal_xi_to_x(double xi, double number_cells, double length_domain, double min_dx) {
    double phase = xi/ number_cells;
    return length_domain * (
        phase - (1.0 /number_cells -  min_dx/ length_domain) *
        std::sin(2.0 * M_PI * phase) /
        std::sin(2.0 * M_PI / number_cells)
    );
}

double sinusoidal_xi_to_dx_dxi(double xi, double number_cells, double length_domain, double min_dx) {
    return length_domain * (1.0/number_cells -  (2.0 * M_PI / number_cells)*(1.0/number_cells - min_dx/ length_domain) *
                    std::cos(2.0 * M_PI * xi / number_cells) / std::sin(2.0 * M_PI / number_cells));
}

non_uniform_domain::non_uniform_domain(int num_cells, double length_domain, int left_boundary_condition, int right_boundary_condition,
    int type, int temp_int_1, int temp_int_2, double temp_double_1, double temp_double_2) {
    this->number_cells = num_cells;
    this->domain_type = 1;
    this->length_domain = length_domain;
    this->left_boundary_condition = left_boundary_condition;
    this->right_boundary_condition = right_boundary_condition;
    this->number_nodes = num_cells + 1;
    this->curvilinear = (temp_int_1 == 1); // curvilinear flag
    // sinusoidal grid generation
    this->grid_nodes.resize(this->number_nodes);
    this->cell_centers.resize(this->number_cells);
    this->dx_dxi.resize(this->number_cells);
    this->grid_nodes[0] = 0.0; // First node
    this->grid_nodes[this->number_cells] = this->length_domain; // Last node
    if (type == 1) {
        // Sinusoidal grid generation
        // Fill the grid using the provided equation
        for (int i = 1; i < this->number_cells; ++i) {
            grid_nodes[i] = sinusoidal_xi_to_x(double(i), double(this->number_cells), this->length_domain, temp_double_1);
        }
        for (int i = 0; i < this->number_cells; ++i) {
            if (!this->curvilinear) {
                this->dx_dxi[i] = this->grid_nodes[i+1] - this->grid_nodes[i];
                this->cell_centers[i] = (this->grid_nodes[i] + this->grid_nodes[i + 1]) * 0.5;
            } else {
                double xi = double(i) + 0.5;
                this->dx_dxi[i] = sinusoidal_xi_to_dx_dxi(xi, double(this->number_cells), this->length_domain, temp_double_1);
                this->cell_centers[i] = sinusoidal_xi_to_x(xi, double(this->number_cells), this->length_domain, temp_double_1); // x at center of logical cell
            }
        }
        this->min_dx = this->dx_dxi[0]; // Minimum cell size
    } else if (type == 2) {
        // half sinusoidal grid generation
        int double_cells = 2 * this->number_cells;
        double double_length = 2.0 * this->length_domain;
        for (int i = 1; i < this->number_cells; ++i) {
            grid_nodes[i] = sinusoidal_xi_to_x(double(i), double_cells, double_length, temp_double_1);
            if (mpi_vars::mpi_rank == 0) {
                std::cout << "grid_nodes[" << i << "] = " << grid_nodes[i] << std::endl;
            }
        }
        for (int i = 0; i < this->number_cells; ++i) {
            if (!this->curvilinear) {
                this->dx_dxi[i] = this->grid_nodes[i+1] - this->grid_nodes[i];
                this->cell_centers[i] = (this->grid_nodes[i] + this->grid_nodes[i + 1]) * 0.5;
            } else {
                double xi = double(i) + 0.5;
                this->dx_dxi[i] = sinusoidal_xi_to_dx_dxi(xi, double_cells, double_length, temp_double_1);
                this->cell_centers[i] = sinusoidal_xi_to_x(xi, double_cells, double_length, temp_double_1);
            }
        }
        this->min_dx = this->dx_dxi[0]; // Minimum cell size
    } else if (type == 3) {
        // sinusoidal with uniform at edges
        int num_uniform_cells = temp_int_2; // Number of uniform cells
        int number_sin_cells = this->number_cells - 2*num_uniform_cells; // Number of sinusoidal cells 
        this->min_dx = temp_double_1 / double(num_uniform_cells); // Uniform cell size
        double length_sin = this->length_domain - 2.0 * temp_double_1; // Length of sinusoidal section
        for (int i = 1; i< num_uniform_cells+1; ++i) {
            this->grid_nodes[i] = i * this->min_dx;
            this->grid_nodes[this->number_cells - i] = this->length_domain - i * this->min_dx;
            this->dx_dxi[i-1] = this->min_dx;
            this->dx_dxi[this->number_cells - i] = this->min_dx;
            this->cell_centers[i-1] = (this->grid_nodes[i-1] + this->grid_nodes[i]) * 0.5;
            this->cell_centers[this->number_cells-i] = (this->grid_nodes[this->number_cells-i+1] + this->grid_nodes[this->number_cells-i]) * 0.5;
        }
        // Find dx_min for sinusoidal function so dx_dxi continous at edges
        double dx_min = (this->min_dx - length_sin / double(number_sin_cells)) / (2.0 * M_PI / double(number_sin_cells) / std::sin(2.0 * M_PI / double(number_sin_cells))) + length_sin / double(number_sin_cells);
        for (int i = 0; i < number_sin_cells; ++i) {
            this->grid_nodes[i + num_uniform_cells + 1] = this->grid_nodes[num_uniform_cells] + sinusoidal_xi_to_x(double(i+1), double(number_sin_cells), length_sin, dx_min);
            if (!this->curvilinear) {
                this->dx_dxi[i + num_uniform_cells] = this->grid_nodes[i + num_uniform_cells+1] - this->grid_nodes[i + num_uniform_cells];
                this->cell_centers[i + num_uniform_cells] = (this->grid_nodes[i + num_uniform_cells+1] + this->grid_nodes[i + num_uniform_cells]) * 0.5;
            } else {
                double xi = double(i) + 0.5;
                this->dx_dxi[i + num_uniform_cells] = sinusoidal_xi_to_dx_dxi(xi, double(number_sin_cells), length_sin, dx_min);
                this->cell_centers[i + num_uniform_cells] = this->grid_nodes[num_uniform_cells] + sinusoidal_xi_to_x(xi, double(number_sin_cells), length_sin, temp_double_1);
            }
        }
    } else if (type == 4) {
        // uniform and half-sinusoidal grid generation
        int num_uniform_cells = temp_int_2; // Number of uniform cells
        int number_sin_cells = this->number_cells - num_uniform_cells; // Number of sinusoidal cells 
        this->min_dx = temp_double_1 / double(num_uniform_cells); // Uniform cell size
        double length_sin = this->length_domain - temp_double_1; // Length of sinusoidal section
        for (int i = 1; i< num_uniform_cells+1; ++i) {
            this->grid_nodes[i] = i * this->min_dx;
            this->dx_dxi[i-1] = this->min_dx;
            this->cell_centers[i-1] = (this->grid_nodes[i-1] + this->grid_nodes[i]) * 0.5;
        }
        // Find dx_min for sinusoidal function so dx_dxi continous at edges
        double dx_min = (this->min_dx - length_sin / double(number_sin_cells)) / (M_PI / double(number_sin_cells) / std::sin(M_PI / double(number_sin_cells))) + length_sin / double(number_sin_cells);
        for (int i = 0; i < number_sin_cells; ++i) {
            this->grid_nodes[i + num_uniform_cells + 1] = this->grid_nodes[num_uniform_cells] + sinusoidal_xi_to_x(double(i+1), 2.0 * double(number_sin_cells), 2.0 * length_sin, dx_min);
            if (!this->curvilinear) {
                this->dx_dxi[i + num_uniform_cells] = this->grid_nodes[i + num_uniform_cells+1] - this->grid_nodes[i + num_uniform_cells];
                this->cell_centers[i + num_uniform_cells] = (this->grid_nodes[i + num_uniform_cells+1] + this->grid_nodes[i + num_uniform_cells]) * 0.5;
            } else {
                double xi = double(i) + 0.5;
                this->dx_dxi[i + num_uniform_cells] = sinusoidal_xi_to_dx_dxi(xi, 2.0 * double(number_sin_cells), 2.0 * length_sin, dx_min);
                this->cell_centers[i + num_uniform_cells] = this->grid_nodes[num_uniform_cells] + sinusoidal_xi_to_x(xi, 2.0 * double(number_sin_cells), 2.0 * length_sin, dx_min);
            }
        }
    }
   
}

void uniform_domain::write_domain(const std::string& filename) {
    // Open file
    if (mpi_vars::mpi_rank == 0) {
        std::ofstream file(filename + "/domain/parameters.dat");
        if (!file) {
            std::cerr << "Error opening file for domain \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Write header (optional)
        file << "type, Number nodes, number cells, left_boundary_condition, right_boundary_condition, length_domain (m) \n";

        file << std::scientific << std::setprecision(8);
        file << this->domain_type << "\t"
        << this->number_nodes << "\t"
        << this->number_cells << "\t"
        << this->left_boundary_condition << "\t"
        << this->right_boundary_condition << "\t"
        << this->length_domain <<
        "\n";

        file.close();

        write_vector_to_binary_file(this->grid_nodes, this->number_nodes, filename + "/domain/grid.dat", 0);
        std::vector<double> temp(1);
        temp[0] = this->min_dx;
        write_vector_to_binary_file(temp, 1, filename + "/domain/dx_dxi.dat", 0);

    }
}

void non_uniform_domain::print_out(){
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "Non-uniform domain with " << this->number_cells << " cells and " << this->number_nodes << " nodes." << std::endl;
        std::cout << "Length of the domain: " << this->length_domain << std::endl;
        std::cout << "Minimum cell size: " << this->min_dx << std::endl;
        std::cout << "Left boundary condition: " << this->left_boundary_condition << std::endl;
        std::cout << "Right boundary condition: " << this->right_boundary_condition << std::endl;
        std::cout << "Curvilinear: " << (this->curvilinear ? "true" : "false") << std::endl;
        std::cout << "--------------------------------" << std::endl;
    }
}

void non_uniform_domain::write_domain(const std::string& filename) {
    // Open file
    if (mpi_vars::mpi_rank == 0) {
        std::ofstream file(filename + "/domain/parameters.dat");
        if (!file) {
            std::cerr << "Error opening file for domain \n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Write header (optional)
        file << "type, Number nodes, number cells, left_boundary_condition, right_boundary_condition, length_domain (m) \n";

        file << std::scientific << std::setprecision(8);
        file << this->domain_type << "\t"
        << this->number_nodes << "\t"
        << this->number_cells << "\t"
        << this->left_boundary_condition << "\t"
        << this->right_boundary_condition << "\t"
        << this->length_domain <<
        "\n";

        file.close();

        write_vector_to_binary_file(this->grid_nodes, this->number_nodes, filename + "/domain/grid.dat", 0);
        write_vector_to_binary_file(this->dx_dxi, this->number_cells, filename + "/domain/dx_dxi.dat", 0);

    }
}

std::unique_ptr<domain> create_domain_from_file(const std::string& filename, int scheme_type) {
    int left_boundary, right_boundary;
    int number_cells;
    double length_domain;
    int type;
    int temp_int_1, temp_int_2;
    double temp_double_1, temp_double_2;
    if (mpi_vars::mpi_rank == 0) {
        std::cout << " "  << std::endl;
        std::cout << "Reading domain inputs: "  << std::endl;
        std::cout << "-------------------------- "  << std::endl;
        std::string line;
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            exit(EXIT_FAILURE);
        }
        
        
        std::getline(file, line);
        std::istringstream iss(line);
        iss >> number_cells;
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> length_domain;
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> type >> temp_int_1 >> temp_int_2 >> temp_double_1 >> temp_double_2;
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> left_boundary >> right_boundary;
        if (left_boundary == 2 && right_boundary == 2) {
            std::cerr << "Error: Both left and right boundary conditions cannot be Neumann (2)." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (left_boundary == 4 && right_boundary == 4) {
            std::cerr << "Error: Both left and right boundary conditions cannot be RF voltage (4)." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (left_boundary == 3 || right_boundary == 3) {
            left_boundary = right_boundary = 3;
        }
        file.close();
    }

    // Broadcast the parameters to all processes
    MPI_Bcast(&type, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&number_cells, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&length_domain, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
    MPI_Bcast(&left_boundary, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&right_boundary, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&temp_double_1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&temp_int_1, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&temp_int_2, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&temp_double_2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (type == 0) {
        // Uniform domain
        return std::make_unique<uniform_domain>(number_cells, length_domain, left_boundary, right_boundary);
    } else { 
        if (scheme_type == 0) {
            std::cout << "Cannot have non-uniform domain with MC-PIC!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        } else {
            // Non-uniform domain               
            return std::make_unique<non_uniform_domain>(number_cells, length_domain, left_boundary, right_boundary,
                type, temp_int_1, temp_int_2, temp_double_1, temp_double_2);
        }
    }
    return nullptr; // This line should never be reached, but added to avoid compiler warnings
}


