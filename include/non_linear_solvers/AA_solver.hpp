#pragma once
#include "non_linear_solver.hpp"
#include "util/math_util.hpp"
#include <functional>
#include "globals/mpi_vars.hpp"
#include <iostream>


class AA_solver : public non_linear_solver { // Anderson accelerated 
public:

    double beta;
    int m_anderson;
    std::vector<double> norm_residual, min_matrix;
    std::vector<std::vector<double>> residual_k, x_k;
    void initialize_diagnostic_files(const std::string& filename) const override; // initialize diagnostic files
    AA_solver(double beta, double eps_r, double eps_a, int max_iterations, int m_anderson, int number_unknowns);
    void print_out() const override;
    void solve(std::vector<double>& x_result, // x_result first with initial guess, then pass actual result
        const std::function<void(std::vector<double>&)>& fixed_point_function) override { //Fixed point function returns next value x_k F(x_k, x_k+1), always uses x_result
        double start_time;
        #pragma omp master
        {   
            start_time = MPI_Wtime();
            this->x_k[0] = x_result;
        }
        #pragma omp barrier
        fixed_point_function(x_result);
        #pragma omp barrier // make sure all threads have finished the fixed point function before proceeding
        #pragma omp master
        {
        
            this->x_k[1] = x_result;
            this->norm_residual[0] = 0.0;
            for (int i = 0; i < this->number_unknowns; i++){
                double diff = this->x_k[1][i] - this->x_k[0][i];
                this->residual_k[0][i] = diff;
                this->norm_residual[0] += diff * diff; // accumulate norm
            }
        
            // if (mpi_vars::mpi_rank == 0) {
            //     std::cout << "Initial residual norm: " << this->norm_residual[0] << std::endl;
            // }
            this->norm_residual[0] = std::sqrt(this->norm_residual[0]);   
        }
        #pragma omp barrier
        double eps_tol = this->eps_r * this->norm_residual[0] + this->eps_a * std::sqrt(double(this->number_unknowns));
        int iter, index, m_k;
        for (iter = 1; iter< this->max_iterations; iter++) {
            if (iter < this->m_anderson) {
                m_k = iter;
            } else {
                m_k = this->m_anderson;
            }
            index = iter % (this->m_anderson+1);
            #pragma omp barrier
            fixed_point_function(x_result); // fixed 
            #pragma omp barrier
            #pragma omp master
            {   
                this->norm_residual[index] = 0.0;
                for (int i = 0; i < this->number_unknowns; i++){
                    double diff = x_result[i] - this->x_k[index][i];
                    this->residual_k[index][i] = diff;
                    this->norm_residual[index] += diff * diff; // accumulate norm
                }
                this->norm_residual[index] = std::sqrt(this->norm_residual[index]);
            } 
            #pragma omp barrier
            
            if (this->norm_residual[index] < eps_tol) {  
                break;
            }
            #pragma omp master
            {
                // if (mpi_vars::mpi_rank == 0) {
                //     std::cout << "Iter is "<< iter << " with norm residual " << this->norm_residual[index] << std::endl;
                //     std::cout << "Index is "<< index << std::endl;
                // }
                for (int j = 0; j < m_k; j++) {
                    int past_indx = (iter - m_k + j) % (this->m_anderson + 1);
                    size_t start_indx = j * this->number_unknowns; // flattened index start
                    for (int i = 0; i < this->number_unknowns; i++) {
                        double diff = this->residual_k[index][i] - this->residual_k[past_indx][i];
                        // if (!std::isfinite(diff)) {
                        //     std::cout << "error" << std::endl;
                        //     MPI_Abort(MPI_COMM_WORLD, 1);
                        // }
                        this->min_matrix[start_indx+i] = diff;
                    }      
                }

                // // solve minimization problem
                // if (mpi_vars::mpi_rank == 0) {
                //     std::cout << "Before alpha " << std::endl;
                // }
                // MPI_Barrier(MPI_COMM_WORLD); // ensure all threads have finished before solving
                std::vector<double> alpha = solveLeastSquaresQR_MKL(min_matrix, residual_k[index], this->number_unknowns, m_k);
                // if (mpi_vars::mpi_rank == 0) {
                //     std::cout << "Went through alpha " << std::endl;
                // }
                
                int next_idx = (index+1) % (this->m_anderson + 1);
                int past_indx = (iter - m_k) % (this->m_anderson + 1);
                double alpha_last = alpha[0];
                // initially set next x_k to first component
                for (int i = 0; i < this->number_unknowns; i++) {
                    this->x_k[next_idx][i] = alpha[0] * (this->beta * residual_k[past_indx][i] + x_k[past_indx][i]);
                }   
                for (int j = 1; j < m_k; j++) {
                    // Sum alphas
                    alpha_last += alpha[j];
                    past_indx = (iter - m_k + j) % (this->m_anderson + 1);
                    for (int i = 0; i < this->number_unknowns; i++) {
                        this->x_k[next_idx][i] += alpha[j] * (this->beta * residual_k[past_indx][i] + x_k[past_indx][i]);
                    }   
                }
                alpha_last = 1.0 - alpha_last; // close coefficients so add to 1
                for (int i = 0; i < this->number_unknowns; i++) {
                    this->x_k[next_idx][i] += alpha_last * (this->beta * residual_k[index][i] + x_k[index][i]); // add current component
                    x_result[i] = this->x_k[next_idx][i];
                }
            }
            #pragma omp barrier

        }
        #pragma omp master
        {   
            double end_time = MPI_Wtime();
            this->solver_time += end_time - start_time;
            this->number_iterations = iter + 1;
            this->accum_residual_norm += this->norm_residual[index];
            this->accum_iterations_count += this->number_iterations;
        }

    }
};


