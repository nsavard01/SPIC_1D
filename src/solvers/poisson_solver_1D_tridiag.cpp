#include "solvers/poisson_solver_1D.hpp"
#include "solvers/poisson_solver_1D_tridiag.hpp"
#include "globals/mpi_vars.hpp"
#include <cmath>

poisson_solver_1D_tridiag::poisson_solver_1D_tridiag(const domain& world) {
    this->number_unknowns = world.number_nodes;
    this->diagonal.resize(this->number_unknowns, 0.0);
    this->upper.resize(this->number_unknowns-1, 0.0);
    this->lower.resize(this->number_unknowns-1, 0.0);
    this->work_space.resize(this->number_unknowns, 0.0);
    if (typeid(world) == typeid(uniform_domain)) {
        double dx = world.min_dx;
        switch (world.left_boundary_condition) {
            case 1:
            case 4:
            case 3:
                this->diagonal[0] = 1.0;
                this->upper[0] = 0.0;
                break;
            case 2:
                this->diagonal[0] = -1.0 / dx;
                this->upper[0] = 1.0 / dx;
                break;
            default:
                throw std::invalid_argument("Invalid left boundary condition.");
        }
    
        switch (world.right_boundary_condition) {
            case 1:
            case 4:
            case 3:
                this->diagonal[this->number_unknowns-1] = 1.0;
                this->lower[this->number_unknowns-2] = 0.0;
                break;
            case 2:
                this->diagonal[this->number_unknowns-1] = -1.0 / dx;
                this->lower[this->number_unknowns-2] = 1.0 / dx;
                break;
            default:
                throw std::invalid_argument("Invalid right boundary condition.");
        }

        for (int i = 1; i < this->number_unknowns - 1; ++i) {
            this->diagonal[i] = -2.0 / dx;
            this->upper[i] = 1.0 / dx;
            this->lower[i-1] = 1.0 / dx;
        }
    } else if (typeid(world) == typeid(non_uniform_domain)) {
        const std::vector<double>& dx_dxi = world.dx_dxi;
        switch (world.left_boundary_condition) {
            case 1:
            case 4:
            case 3:
                this->diagonal[0] = 1.0;
                this->upper[0] = 0.0;
                break;
            case 2:
                this->diagonal[0] = -1.0 / dx_dxi[0];
                this->upper[0] = 1.0 / dx_dxi[0];
                break;
            default:
                throw std::invalid_argument("Invalid left boundary condition.");
        }
    
        switch (world.right_boundary_condition) {
            case 1:
            case 4:
            case 3:
                this->diagonal[this->number_unknowns-1] = 1.0;
                this->lower[this->number_unknowns-2] = 0.0;
                break;
            case 2:
                this->diagonal[this->number_unknowns-1] = -1.0 / dx_dxi[this->number_unknowns-2];
                this->lower[this->number_unknowns-2] = 1.0 / dx_dxi[this->number_unknowns-2];
                break;
            default:
                throw std::invalid_argument("Invalid right boundary condition.");
        }

        for (int i = 1; i < this->number_unknowns - 1; ++i) {
            this->diagonal[i] = -(1.0 / dx_dxi[i-1] + 1.0 / dx_dxi[i]);
            this->upper[i] = 1.0 / dx_dxi[i];
            this->lower[i-1] = 1.0 / dx_dxi[i-1];
        }
    }

    // if (mpi_vars::mpi_rank == 0) {
    //     for (int i = 0; i < this->number_unknowns; ++i) {
    //         std::cout << "Diagonal[" << i << "] = " << this->diagonal[i] << std::endl;
    //         if (i < this->number_unknowns - 1) {
    //             std::cout << "Upper[" << i << "] = " << this->upper[i] << std::endl;
    //             std::cout << "Lower[" << i << "] = " << this->lower[i] << std::endl;
    //         }
    //     }
    // }
    
    
}

void poisson_solver_1D_tridiag::solve(std::vector<double>& solution, std::vector<double>& source_term) {
    // Solve the Poisson equation using the Thomas algorithm
    // Note one can have same reference for solution and source_term if you want to overwrite the source_term
    double m;
    this->work_space[0] = this->upper[0] / this->diagonal[0];
    solution[0] = source_term[0] / this->diagonal[0];
    for (int i = 1; i < this->number_unknowns-1; i++) {
        m = this->diagonal[i] - this->lower[i-1] * this->work_space[i-1];
        this->work_space[i] = this->upper[i] / m;
        solution[i] = (source_term[i] - this->lower[i-1] * solution[i-1]) / m;
    }

    m = this->diagonal[this->number_unknowns-1] - this->lower[this->number_unknowns-2] * this->work_space[this->number_unknowns-2];
    solution[this->number_unknowns-1] = (source_term[this->number_unknowns-1] - this->lower[this->number_unknowns-2] * solution[this->number_unknowns-2]) / m; 
    for (int i = this->number_unknowns-2; i >= 0; i--) {
        solution[i] = solution[i] - this->work_space[i] * solution[i+1];
    }
}

double poisson_solver_1D_tridiag::norm_error(const std::vector<double>& solution, const std::vector<double>& source_term) const {
    // Solve the Poisson equation using the Thomas algorithm
    // Note one can have same reference for solution and source_term if you want to overwrite the source_term
    double error = 0.0;
    double res;
    if (source_term[0] != 0.0) {
        res = (this->diagonal[0] * solution[0] + this->upper[0] * solution[1])/ source_term[0] - 1.0; 
        error += res * res;
    }
   
    for (int i = 1; i < this->number_unknowns-1; i++) {
        if (source_term[i] != 0.0) {
            res = (this->lower[i-1] * solution[i-1] + this->diagonal[i] * solution[i] + this->upper[i] * solution[i+1]) / source_term[i] - 1.0;
            error += res * res;
        }
    }
    
    if (source_term[this->number_unknowns-1] != 0.0) {
        res = (this->lower[this->number_unknowns-2] * solution[this->number_unknowns-2] + this->diagonal[this->number_unknowns-1] * solution[this->number_unknowns-1]) / source_term[this->number_unknowns-1] - 1.0;
        error += res * res;
    }
    
    
    return std::sqrt(error / static_cast<double>(this->number_unknowns));   
}


