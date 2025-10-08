
#include <vector>
#include "particles/charged_particle.hpp"
#include "globals/mpi_vars.hpp"
#include "globals/constants.hpp"
#include "globals/write_functions.hpp"
#include <omp.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <mpi.h>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <dirent.h>
#include "rand_gen/maxwell_generator.hpp"




void charged_particle::ES_push_deposit_INGP_uniform(const int thread_id, double del_t, const std::vector<double>& E_field, std::vector<double>& work_space, 
    const double inv_dx, const int left_boundary, const int right_boundary, const int number_cells) {
   
    size_t last_idx = this->number_particles[thread_id][0];
    const std::vector<double>& xi_local = this->xi[thread_id];
    const std::vector<double>& v_x_local = this->v_x[thread_id];
    double v_x_i, v_x_f, xi_i, xi_f;
    double E_field_local; // E_x in SI units (physical space)
    double del_tau, time_passed, accel, d_interp;
    int xi_cell;
    bool in_cell_bool, equal_v_sign_bool, future_boundary_bool, del_part;
    int v_sign;
    int xi_boundary;
    const double t_tol = del_t * 1e-10;
    const double q_over_m = this->q_over_m;
    for (size_t part_indx= 0; part_indx < last_idx; part_indx++){

        xi_i = xi_local[part_indx];
        v_x_i = v_x_local[part_indx];
        
        v_sign = (v_x_i > 0) - (v_x_i < 0);
        xi_cell = int(xi_i + v_sign*1e-12); // in case lands on exact boundary, has happened before...
        del_tau = del_t;
        time_passed = 0.0;
        del_part = false;
        while (del_tau > t_tol) {
            E_field_local = E_field[xi_cell];
            accel = q_over_m * E_field_local;
            v_x_f = v_x_i + accel * del_tau;
            xi_f = xi_i + 0.5 * (v_x_i + v_x_f) * del_tau * inv_dx;
        
            // check if particle outside cell or flipped direction
            in_cell_bool = (int(xi_f+1) == xi_cell+1);
            equal_v_sign_bool = ((v_x_i > 0) == (v_x_f > 0));
            future_boundary_bool = (!in_cell_bool) || (!equal_v_sign_bool);

            if (future_boundary_bool) {
                
                if (v_x_i > 0) {
                    xi_boundary = xi_cell + 1;
                } else {
                    xi_boundary = xi_cell;
                }
                double diff_PE = 2.0 * accel * (double(xi_boundary) - xi_i) / inv_dx;
                double v_i_sqr = v_x_i*v_x_i;
                // will particle reach boundary?
                future_boundary_bool = (diff_PE > -v_i_sqr);
                if (future_boundary_bool) {
                    xi_f = double(xi_boundary);
                    v_x_f = v_sign * std::sqrt(diff_PE + v_i_sqr);
                    // select method to minimize error in del_tau
                    if (std::abs(v_x_f - v_x_i) > std::abs(v_x_f + v_x_i)) {
                        del_tau = (v_x_f - v_x_i)/accel;
                    } else {
                        del_tau = 2.0 * (xi_f - xi_i)/ (v_x_f + v_x_i) / inv_dx;
                    }
                } else if (!in_cell_bool) {
                    future_boundary_bool = true;
                    xi_boundary = 2 * xi_cell + 1 - xi_boundary;
                    xi_f = double(xi_boundary);
                    diff_PE = 2.0 * accel * (xi_f - xi_i) / inv_dx;
                    v_x_f = - v_sign * std::sqrt(diff_PE + v_i_sqr);
                    del_tau = (v_x_f - v_x_i)/accel;
                }

            }
            
            if (future_boundary_bool) {
                v_sign = (v_x_f > 0) - (v_x_f < 0);
                xi_cell = xi_cell + v_sign;
                if (xi_boundary == 0) {
                    switch (left_boundary){
                        case 1:
                        case 4:
                            xi_cell = 0;
                            del_part = true;
                            break;
                        case 2:
                            xi_cell = 0;
                            v_x_f = - v_x_f;
                            v_sign = -v_sign;
                            break;
                        case 3:
                            xi_f = double(number_cells);
                            xi_cell = number_cells-1;
                            break;
                    }
                } else if (xi_boundary == number_cells) {
                    switch (right_boundary){
                        case 1:
                        case 4:
                            xi_cell = number_cells-1;
                            del_part = true;
                            break;
                        case 2:
                            xi_cell = number_cells-1;
                            v_x_f = -v_x_f;
                            v_sign = -v_sign;
                            break;
                        case 3:
                            xi_f = 0.0;
                            xi_cell = 0;
                            break;
                    }
    
                }
                if (del_part) {
                    // particle is lost, no need to continue
                    break;
                }
            }
            
            time_passed += del_tau;
            del_tau = del_t - time_passed;
            xi_i = xi_f;
            v_x_i = v_x_f;
            
        }
        d_interp = xi_f - xi_cell;
        work_space[xi_cell] += (1.0 - d_interp);
        work_space[xi_cell+1] += d_interp;
        
    }
    for (int inj_indx = 0; inj_indx < this->number_unique_injections; inj_indx++) {
        size_t start_indx = last_idx;
        const size_t number_inject_local = this->number_particles_injected[thread_id][inj_indx];
        const std::vector<double>& del_t_array = this->time_step_injected[thread_id][inj_indx];
        last_idx = start_indx + number_inject_local;
        double del_t_local;
        for (size_t part_indx = start_indx; part_indx < last_idx; part_indx++){
            xi_i = xi_local[part_indx];
            v_x_i = v_x_local[part_indx];
            del_t_local = del_t_array[part_indx-start_indx];

            v_sign = (v_x_i > 0) - (v_x_i < 0);
            xi_cell = int(xi_i + v_sign*1e-12); // in case lands on exact boundary, has happened before...
            del_tau = del_t_local;
            time_passed = 0.0;
            del_part = false;
            while (del_tau > t_tol) {
                E_field_local = E_field[xi_cell];
                accel = q_over_m * E_field_local;
                v_x_f = v_x_i + accel * del_tau;
                xi_f = xi_i + 0.5 * (v_x_i + v_x_f) * del_tau * inv_dx;
            
                // check if particle outside cell or flipped direction
                in_cell_bool = (int(xi_f+1) == xi_cell+1);
                equal_v_sign_bool = ((v_x_i > 0) == (v_x_f > 0));
                future_boundary_bool = (!in_cell_bool) || (!equal_v_sign_bool);

                if (future_boundary_bool) {
                    
                    if (v_x_i > 0) {
                        xi_boundary = xi_cell + 1;
                    } else {
                        xi_boundary = xi_cell;
                    }
                    double diff_PE = 2.0 * accel * (double(xi_boundary) - xi_i) / inv_dx;
                    double v_i_sqr = v_x_i*v_x_i;
                    // will particle reach boundary?
                    future_boundary_bool = (diff_PE > -v_i_sqr);
                    if (future_boundary_bool) {
                        xi_f = double(xi_boundary);
                        v_x_f = v_sign * std::sqrt(diff_PE + v_i_sqr);
                        // select method to minimize error in del_tau
                        if (std::abs(v_x_f - v_x_i) > std::abs(v_x_f + v_x_i)) {
                            del_tau = (v_x_f - v_x_i)/accel;
                        } else {
                            del_tau = 2.0 * (xi_f - xi_i)/ (v_x_f + v_x_i) / inv_dx;
                        }
                    } else if (!in_cell_bool) {
                        future_boundary_bool = true;
                        xi_boundary = 2 * xi_cell + 1 - xi_boundary;
                        xi_f = double(xi_boundary);
                        diff_PE = 2.0 * accel * (xi_f - xi_i) / inv_dx;
                        v_x_f = - v_sign * std::sqrt(diff_PE + v_i_sqr);
                        del_tau = (v_x_f - v_x_i)/accel;
                    }

                }
                
                if (future_boundary_bool) {
                    v_sign = (v_x_f > 0) - (v_x_f < 0);
                    xi_cell = xi_cell + v_sign;
                    if (xi_boundary == 0) {
                        switch (left_boundary){
                            case 1:
                            case 4:
                                xi_cell = 0;
                                del_part = true;
                                break;
                            case 2:
                                xi_cell = 0;
                                v_x_f = - v_x_f;
                                v_sign = -v_sign;
                                break;
                            case 3:
                                xi_f = double(number_cells);
                                xi_cell = number_cells-1;
                                break;
                        }
                    } else if (xi_boundary == number_cells) {
                        switch (right_boundary){
                            case 1:
                            case 4:
                                xi_cell = number_cells-1;
                                del_part = true;
                                break;
                            case 2:
                                xi_cell = number_cells-1;
                                v_x_f = -v_x_f;
                                v_sign = -v_sign;
                                break;
                            case 3:
                                xi_f = 0.0;
                                xi_cell = 0;
                                break;
                        }
        
                    }
                    if (del_part) {
                        // particle is lost, no need to continue
                        break;
                    }
                }
                
                time_passed += del_tau;
                del_tau = del_t_local - time_passed;
                xi_i = xi_f;
                v_x_i = v_x_f;
                
            }
            d_interp = xi_f - xi_cell;
            work_space[xi_cell] += (1.0 - d_interp);
            work_space[xi_cell+1] += d_interp;
        }
    }
}

void charged_particle::ES_push_INGP_uniform(const int thread_id, double del_t, const std::vector<double>& E_field, int& number_sub_steps, 
    const double inv_dx, const int left_boundary, const int right_boundary, const int number_cells) {
   
    size_t last_idx = this->number_particles[thread_id][0];
    std::vector<double>& xi_local = this->xi[thread_id];
    std::vector<double>& v_x_local = this->v_x[thread_id];
    bool use_vy = (this->number_velocity_coordinates > 1);
    bool use_vz = (this->number_velocity_coordinates > 2);
    bool use_y = (this->number_space_coordinates > 1);
    bool use_z = (this->number_space_coordinates > 2);
    double v_x_i, v_x_f, xi_i, xi_f, v_y = 0.0, v_z = 0.0;
    double E_field_local; // E_x in SI units (physical space)
    double del_tau, time_passed, accel;
    int xi_cell;
    bool in_cell_bool, equal_v_sign_bool, future_boundary_bool, del_part;
    int v_sign;
    int xi_boundary;
    const double t_tol = del_t * 1e-10;
    const double q_over_m = this->q_over_m;
    size_t space_delete = 0;
    for (size_t part_indx= 0; part_indx < last_idx; part_indx++){

        xi_i = xi_local[part_indx];
        v_x_i = v_x_local[part_indx];
        
        v_sign = (v_x_i > 0) - (v_x_i < 0);
        xi_cell = int(xi_i + v_sign*1e-12); // in case lands on exact boundary, has happened before...
        del_tau = del_t;
        time_passed = 0.0;
        del_part = false;
        if (use_vy) {
            v_y = this->v_y[thread_id][part_indx];
        }
        if (use_vz){
            v_z = this->v_z[thread_id][part_indx];
        }
        while (del_tau > t_tol) {
            number_sub_steps++;
            E_field_local = E_field[xi_cell];
            accel = q_over_m * E_field_local;
            v_x_f = v_x_i + accel * del_tau;
            xi_f = xi_i + 0.5 * (v_x_i + v_x_f) * del_tau * inv_dx;
            
            // check if particle outside cell or flipped direction
            in_cell_bool = (int(xi_f+1) == xi_cell+1);
            equal_v_sign_bool = ((v_x_i > 0) == (v_x_f > 0));
            future_boundary_bool = (!in_cell_bool) || (!equal_v_sign_bool);

            if (future_boundary_bool) {
                
                if (v_x_i > 0) {
                    xi_boundary = xi_cell + 1;
                } else {
                    xi_boundary = xi_cell;
                }
                double diff_PE = 2.0 * accel * (double(xi_boundary) - xi_i) / inv_dx;
                double v_i_sqr = v_x_i*v_x_i;
                // will particle reach boundary?
                future_boundary_bool = (diff_PE > -v_i_sqr);
                if (future_boundary_bool) {
                    xi_f = double(xi_boundary);
                    v_x_f = v_sign * std::sqrt(diff_PE + v_i_sqr);
                    // select method to minimize error in del_tau
                    if (std::abs(v_x_f - v_x_i) > std::abs(v_x_f + v_x_i)) {
                        del_tau = (v_x_f - v_x_i)/accel;
                    } else {
                        del_tau = 2.0 * (xi_f - xi_i)/ (v_x_f + v_x_i) / inv_dx;
                    }
                } else if (!in_cell_bool) {
                    future_boundary_bool = true;
                    xi_boundary = 2 * xi_cell + 1 - xi_boundary;
                    xi_f = double(xi_boundary);
                    diff_PE = 2.0 * accel * (xi_f - xi_i) / inv_dx;
                    v_x_f = - v_sign * std::sqrt(diff_PE + v_i_sqr);
                    del_tau = (v_x_f - v_x_i)/accel;
                }

            }

            if (future_boundary_bool) {
                v_sign = (v_x_f > 0) - (v_x_f < 0);
                xi_cell = xi_cell + v_sign;
                if (xi_boundary == 0) {
                    switch (left_boundary){
                        case 1:
                        case 4:
                            this->energy_loss[thread_id][0] += v_x_f*v_x_f + v_y*v_y + v_z*v_z;
                            this->wall_loss[thread_id][0]++;
                            this->momentum_loss[thread_id][0][0] += v_x_f;
                            this->momentum_loss[thread_id][0][1] += v_y;
                            this->momentum_loss[thread_id][0][2] += v_z;
                            del_part = true;
                            break;
                        case 2:
                            xi_cell = 0;
                            v_x_f = - v_x_f;
                            v_sign = -v_sign;
                            break;
                        case 3:
                            xi_f = double(number_cells);
                            xi_cell = number_cells-1;
                            break;
                    }
                } else if (xi_boundary == number_cells) {
                    switch (right_boundary){
                        case 1:
                        case 4:
                            this->energy_loss[thread_id][1] += v_x_f*v_x_f + v_y*v_y + v_z*v_z;
                            this->wall_loss[thread_id][1]++;
                            this->momentum_loss[thread_id][1][0] += v_x_f;
                            this->momentum_loss[thread_id][1][1] += v_y;
                            this->momentum_loss[thread_id][1][2] += v_z;
                            del_part = true;
                            break;
                        case 2:
                            xi_cell = number_cells-1;
                            v_x_f = -v_x_f;
                            v_sign = -v_sign;
                            break;
                        case 3:
                            xi_f = 0.0;
                            xi_cell = 0;
                            break;
                    }
    
                }
                if (del_part) {
                    // particle is lost, no need to continue
                    break;
                }
            }
            
            time_passed += del_tau;
            del_tau = del_t - time_passed;
            xi_i = xi_f;
            v_x_i = v_x_f;
            
        }

        if (!del_part) {
            size_t new_idx = part_indx-space_delete;
            xi_local[new_idx] = xi_f;
            v_x_local[new_idx] = v_x_f;
            
            if (use_vy) {
                this->v_y[thread_id][new_idx] = v_y;
            }
            if (use_vz){
                this->v_z[thread_id][new_idx] = v_z;
            }
            if (use_y) {
                this->y[thread_id][new_idx] = this->y[thread_id][part_indx];
            }
            if (use_z) {
                this->z[thread_id][new_idx] = this->z[thread_id][part_indx];
            }
        } else {
            space_delete++;
        }
        
    }
    for (int inj_indx = 0; inj_indx < this->number_unique_injections; inj_indx++) {
        size_t start_indx = last_idx;
        const size_t number_inject_local = this->number_particles_injected[thread_id][inj_indx];
        const std::vector<double>& del_t_array = this->time_step_injected[thread_id][inj_indx];
        last_idx = start_indx + number_inject_local;
        double del_t_local;
        for (size_t part_indx = start_indx; part_indx < last_idx; part_indx++){
            xi_i = xi_local[part_indx];
            v_x_i = v_x_local[part_indx];
            del_t_local = del_t_array[part_indx-start_indx];

            v_sign = (v_x_i > 0) - (v_x_i < 0);
            xi_cell = int(xi_i + v_sign*1e-12); // in case lands on exact boundary, has happened before...
            del_tau = del_t_local;
            time_passed = 0.0;
            del_part = false;
            if (use_vy) {
                v_y = this->v_y[thread_id][part_indx];
            }
            if (use_vz){
                v_z = this->v_z[thread_id][part_indx];
            }
            while (del_tau > t_tol) {
                number_sub_steps++;
                E_field_local = E_field[xi_cell];
                accel = q_over_m * E_field_local;
                v_x_f = v_x_i + accel * del_tau;
                xi_f = xi_i + 0.5 * (v_x_i + v_x_f) * del_tau * inv_dx;
                
                // check if particle outside cell or flipped direction
                in_cell_bool = (int(xi_f+1) == xi_cell+1);
                equal_v_sign_bool = ((v_x_i > 0) == (v_x_f > 0));
                future_boundary_bool = (!in_cell_bool) || (!equal_v_sign_bool);

                if (future_boundary_bool) {
                    
                    if (v_x_i > 0) {
                        xi_boundary = xi_cell + 1;
                    } else {
                        xi_boundary = xi_cell;
                    }
                    double diff_PE = 2.0 * accel * (double(xi_boundary) - xi_i) / inv_dx;
                    double v_i_sqr = v_x_i*v_x_i;
                    // will particle reach boundary?
                    future_boundary_bool = (diff_PE > -v_i_sqr);
                    if (future_boundary_bool) {
                        xi_f = double(xi_boundary);
                        v_x_f = v_sign * std::sqrt(diff_PE + v_i_sqr);
                        // select method to minimize error in del_tau
                        if (std::abs(v_x_f - v_x_i) > std::abs(v_x_f + v_x_i)) {
                            del_tau = (v_x_f - v_x_i)/accel;
                        } else {
                            del_tau = 2.0 * (xi_f - xi_i)/ (v_x_f + v_x_i) / inv_dx;
                        }
                    } else if (!in_cell_bool) {
                        future_boundary_bool = true;
                        xi_boundary = 2 * xi_cell + 1 - xi_boundary;
                        xi_f = double(xi_boundary);
                        diff_PE = 2.0 * accel * (xi_f - xi_i) / inv_dx;
                        v_x_f = - v_sign * std::sqrt(diff_PE + v_i_sqr);
                        del_tau = (v_x_f - v_x_i)/accel;
                    }

                }

                if (future_boundary_bool) {
                    v_sign = (v_x_f > 0) - (v_x_f < 0);
                    xi_cell = xi_cell + v_sign;
                    if (xi_boundary == 0) {
                        switch (left_boundary){
                            case 1:
                            case 4:
                                this->energy_loss[thread_id][0] += v_x_f*v_x_f + v_y*v_y + v_z*v_z;
                                this->wall_loss[thread_id][0]++;
                                this->momentum_loss[thread_id][0][0] += v_x_f;
                                this->momentum_loss[thread_id][0][1] += v_y;
                                this->momentum_loss[thread_id][0][2] += v_z;
                                del_part = true;
                                break;
                            case 2:
                                xi_cell = 0;
                                v_x_f = - v_x_f;
                                v_sign = -v_sign;
                                break;
                            case 3:
                                xi_f = double(number_cells);
                                xi_cell = number_cells-1;
                                break;
                        }
                    } else if (xi_boundary == number_cells) {
                        switch (right_boundary){
                            case 1:
                            case 4:
                                this->energy_loss[thread_id][1] += v_x_f*v_x_f + v_y*v_y + v_z*v_z;
                                this->wall_loss[thread_id][1]++;
                                this->momentum_loss[thread_id][1][0] += v_x_f;
                                this->momentum_loss[thread_id][1][1] += v_y;
                                this->momentum_loss[thread_id][1][2] += v_z;
                                del_part = true;
                                break;
                            case 2:
                                xi_cell = number_cells-1;
                                v_x_f = -v_x_f;
                                v_sign = -v_sign;
                                break;
                            case 3:
                                xi_f = 0.0;
                                xi_cell = 0;
                                break;
                        }
        
                    }
                    if (del_part) {
                        // particle is lost, no need to continue
                        break;
                    }
                }
                
                time_passed += del_tau;
                del_tau = del_t_local - time_passed;
                xi_i = xi_f;
                v_x_i = v_x_f;
                
            }

            if (!del_part) {
                size_t new_idx = part_indx-space_delete;
                xi_local[new_idx] = xi_f;
                v_x_local[new_idx] = v_x_f;
                
                if (use_vy) {
                    this->v_y[thread_id][new_idx] = v_y;
                }
                if (use_vz){
                    this->v_z[thread_id][new_idx] = v_z;
                }
                if (use_y) {
                    this->y[thread_id][new_idx] = this->y[thread_id][part_indx];
                }
                if (use_z) {
                    this->z[thread_id][new_idx] = this->z[thread_id][part_indx];
                }
            } else {
                space_delete++;
            }
        }
    }
    this->number_particles[thread_id][0] = (last_idx - space_delete);
    this->number_collidable_particles[thread_id][0] = (last_idx - space_delete);
}

void charged_particle::ES_push_deposit_INGP_non_uniform(const int thread_id, double del_t, const std::vector<double>& E_field, std::vector<double>& work_space, 
    const std::vector<double>& dx_dxi, const int left_boundary, const int right_boundary, const int number_cells) {
   
    size_t last_idx = this->number_particles[thread_id][0];
    const std::vector<double>& xi_local = this->xi[thread_id];
    const std::vector<double>& v_x_local = this->v_x[thread_id];
    double v_x_i, v_x_f, xi_i, xi_f;
    double E_field_local; // E_x in SI units (physical space)
    double del_tau, time_passed, accel, d_interp;
    int xi_cell;
    bool in_cell_bool, equal_v_sign_bool, future_boundary_bool, del_part;
    int v_sign;
    int xi_boundary;
    const double t_tol = del_t * 1e-10;
    double dx;
    const double q_over_m = this->q_over_m;
    for (size_t part_indx= 0; part_indx < last_idx; part_indx++){

        xi_i = xi_local[part_indx];
        v_x_i = v_x_local[part_indx];
        v_sign = (v_x_i > 0) - (v_x_i < 0);
        xi_cell = int(xi_i + v_sign*1e-12); // in case lands on exact boundary, has happened before...
        del_tau = del_t;
        time_passed = 0.0;
        del_part = false;
        while (del_tau > t_tol) {
            E_field_local = E_field[xi_cell];
            dx = dx_dxi[xi_cell];
            accel = q_over_m * E_field_local;
            v_x_f = v_x_i + accel * del_tau;
            xi_f = xi_i + 0.5 * (v_x_i + v_x_f) * del_tau / dx;
            
            // check if particle outside cell or flipped direction
            in_cell_bool = (int(xi_f+1) == xi_cell+1);
            equal_v_sign_bool = ((v_x_i > 0) == (v_x_f > 0));
            future_boundary_bool = (!in_cell_bool) || (!equal_v_sign_bool);

            if (future_boundary_bool) {
                
                if (v_x_i > 0) {
                    xi_boundary = xi_cell + 1;
                } else {
                    xi_boundary = xi_cell;
                }
                double diff_PE = 2.0 * accel * (double(xi_boundary) - xi_i) * dx;
                double v_i_sqr = v_x_i*v_x_i;
                // will particle reach boundary?
                future_boundary_bool = (diff_PE > -v_i_sqr);
                if (future_boundary_bool) {
                    xi_f = double(xi_boundary);
                    v_x_f = v_sign * std::sqrt(diff_PE + v_i_sqr);
                    // select method to minimize error in del_tau
                    if (std::abs(v_x_f - v_x_i) > std::abs(v_x_f + v_x_i)) {
                        del_tau = (v_x_f - v_x_i)/accel;
                    } else {
                        del_tau = 2.0 * (xi_f - xi_i)/ (v_x_f + v_x_i) * dx;
                    }
                } else if (!in_cell_bool) {
                    future_boundary_bool = true;
                    xi_boundary = 2 * xi_cell + 1 - xi_boundary;
                    xi_f = double(xi_boundary);
                    diff_PE = 2.0 * accel * (xi_f - xi_i) * dx;
                    v_x_f = - v_sign * std::sqrt(diff_PE + v_i_sqr);
                    del_tau = (v_x_f - v_x_i)/accel;
                }

            }

            if (future_boundary_bool) {
                v_sign = (v_x_f > 0) - (v_x_f < 0);
                xi_cell = xi_cell + v_sign;
                if (xi_boundary == 0) {
                    switch (left_boundary){
                        case 1:
                        case 4:
                            xi_cell = 0;
                            del_part = true;
                            break;
                        case 2:
                            xi_cell = 0;
                            v_x_f = - v_x_f;
                            v_sign = -v_sign;
                            break;
                        case 3:
                            xi_f = double(number_cells);
                            xi_cell = number_cells-1;
                            break;
                    }
                } else if (xi_boundary == number_cells) {
                    switch (right_boundary){
                        case 1:
                        case 4:
                            xi_cell = number_cells-1;
                            del_part = true;
                            break;
                        case 2:
                            xi_cell = number_cells-1;
                            v_x_f = -v_x_f;
                            v_sign = -v_sign;
                            break;
                        case 3:
                            xi_f = 0.0;
                            xi_cell = 0;
                            break;
                    }
    
                }
                if (del_part) {
                    // particle is lost, no need to continue
                    break;
                }
            }
            
            time_passed += del_tau;
            del_tau = del_t - time_passed;
            xi_i = xi_f;
            v_x_i = v_x_f;
            
        }
        d_interp = xi_f - xi_cell;
        work_space[xi_cell] += (1.0 - d_interp);
        work_space[xi_cell+1] += d_interp;
        
    }
    for (int inj_indx = 0; inj_indx < this->number_unique_injections; inj_indx++) {
        size_t start_indx = last_idx;
        const size_t number_inject_local = this->number_particles_injected[thread_id][inj_indx];
        const std::vector<double>& del_t_array = this->time_step_injected[thread_id][inj_indx];
        last_idx = start_indx + number_inject_local;
        double del_t_local;
        for (size_t part_indx = start_indx; part_indx < last_idx; part_indx++){
            xi_i = xi_local[part_indx];
            v_x_i = v_x_local[part_indx];
            del_t_local = del_t_array[part_indx-start_indx];
            v_sign = (v_x_i > 0) - (v_x_i < 0);
            xi_cell = int(xi_i + v_sign*1e-12); // in case lands on exact boundary, has happened before...
            del_tau = del_t_local;
            time_passed = 0.0;
            del_part = false;
            while (del_tau > t_tol) {
                E_field_local = E_field[xi_cell];
                dx = dx_dxi[xi_cell];
                accel = q_over_m * E_field_local;
                v_x_f = v_x_i + accel * del_tau;
                xi_f = xi_i + 0.5 * (v_x_i + v_x_f) * del_tau / dx;
                
                // check if particle outside cell or flipped direction
                in_cell_bool = (int(xi_f+1) == xi_cell+1);
                equal_v_sign_bool = ((v_x_i > 0) == (v_x_f > 0));
                future_boundary_bool = (!in_cell_bool) || (!equal_v_sign_bool);

                if (future_boundary_bool) {
                    
                    if (v_x_i > 0) {
                        xi_boundary = xi_cell + 1;
                    } else {
                        xi_boundary = xi_cell;
                    }
                    double diff_PE = 2.0 * accel * (double(xi_boundary) - xi_i) * dx;
                    double v_i_sqr = v_x_i*v_x_i;
                    // will particle reach boundary?
                    future_boundary_bool = (diff_PE > -v_i_sqr);
                    if (future_boundary_bool) {
                        xi_f = double(xi_boundary);
                        v_x_f = v_sign * std::sqrt(diff_PE + v_i_sqr);
                        // select method to minimize error in del_tau
                        if (std::abs(v_x_f - v_x_i) > std::abs(v_x_f + v_x_i)) {
                            del_tau = (v_x_f - v_x_i)/accel;
                        } else {
                            del_tau = 2.0 * (xi_f - xi_i)/ (v_x_f + v_x_i) * dx;
                        }
                    } else if (!in_cell_bool) {
                        future_boundary_bool = true;
                        xi_boundary = 2 * xi_cell + 1 - xi_boundary;
                        xi_f = double(xi_boundary);
                        diff_PE = 2.0 * accel * (xi_f - xi_i) * dx;
                        v_x_f = - v_sign * std::sqrt(diff_PE + v_i_sqr);
                        del_tau = (v_x_f - v_x_i)/accel;
                    }

                }

                if (future_boundary_bool) {
                    v_sign = (v_x_f > 0) - (v_x_f < 0);
                    xi_cell = xi_cell + v_sign;
                    if (xi_boundary == 0) {
                        switch (left_boundary){
                            case 1:
                            case 4:
                                xi_cell = 0;
                                del_part = true;
                                break;
                            case 2:
                                xi_cell = 0;
                                v_x_f = - v_x_f;
                                v_sign = -v_sign;
                                break;
                            case 3:
                                xi_f = double(number_cells);
                                xi_cell = number_cells-1;
                                break;
                        }
                    } else if (xi_boundary == number_cells) {
                        switch (right_boundary){
                            case 1:
                            case 4:
                                xi_cell = number_cells-1;
                                del_part = true;
                                break;
                            case 2:
                                xi_cell = number_cells-1;
                                v_x_f = -v_x_f;
                                v_sign = -v_sign;
                                break;
                            case 3:
                                xi_f = 0.0;
                                xi_cell = 0;
                                break;
                        }
        
                    }
                    if (del_part) {
                        // particle is lost, no need to continue
                        break;
                    }
                }
                
                time_passed += del_tau;
                del_tau = del_t_local - time_passed;
                xi_i = xi_f;
                v_x_i = v_x_f;
                
            }
            d_interp = xi_f - xi_cell;
            work_space[xi_cell] += (1.0 - d_interp);
            work_space[xi_cell+1] += d_interp;
        }
    }
}


void charged_particle::ES_push_INGP_non_uniform(const int thread_id, double del_t, const std::vector<double>& E_field,
    int& number_sub_steps, const std::vector<double>& dx_dxi, const int left_boundary, const int right_boundary, const int number_cells) {
   
    size_t last_idx = this->number_particles[thread_id][0];
    std::vector<double>& xi_local = this->xi[thread_id];
    std::vector<double>& v_x_local = this->v_x[thread_id];
    bool use_vy = (this->number_velocity_coordinates > 1);
    bool use_vz = (this->number_velocity_coordinates > 2);
    bool use_y = (this->number_space_coordinates > 1);
    bool use_z = (this->number_space_coordinates > 2);
    double v_x_i, v_x_f, xi_i, xi_f, v_y = 0.0, v_z = 0.0;
    double E_field_local; // E_x in SI units (physical space)
    double del_tau, time_passed, accel;
    int xi_cell;
    bool in_cell_bool, equal_v_sign_bool, future_boundary_bool, del_part;
    int v_sign;
    int xi_boundary;
    const double t_tol = del_t * 1e-10;
    double dx;
    const double q_over_m = this->q_over_m;
    size_t space_delete = 0;
    for (size_t part_indx= 0; part_indx < last_idx; part_indx++){

        xi_i = xi_local[part_indx];
        v_x_i = v_x_local[part_indx];
        v_sign = (v_x_i > 0) - (v_x_i < 0);
        xi_cell = int(xi_i + v_sign*1e-12); // in case lands on exact boundary, has happened before...
        del_tau = del_t;
        time_passed = 0.0;
        del_part = false;
        if (use_vy) {
            v_y = this->v_y[thread_id][part_indx];
        }
        if (use_vz){
            v_z = this->v_z[thread_id][part_indx];
        }
        while (del_tau > t_tol) {
            number_sub_steps++;
            E_field_local = E_field[xi_cell];
            dx = dx_dxi[xi_cell];
            accel = q_over_m * E_field_local;
            v_x_f = v_x_i + accel * del_tau;
            xi_f = xi_i + 0.5 * (v_x_i + v_x_f) * del_tau / dx;
            
            // check if particle outside cell or flipped direction
            in_cell_bool = (int(xi_f+1) == xi_cell+1);
            equal_v_sign_bool = ((v_x_i > 0) == (v_x_f > 0));
            future_boundary_bool = (!in_cell_bool) || (!equal_v_sign_bool);

            if (future_boundary_bool) {
                
                if (v_x_i > 0) {
                    xi_boundary = xi_cell + 1;
                } else {
                    xi_boundary = xi_cell;
                }
                double diff_PE = 2.0 * accel * (double(xi_boundary) - xi_i) * dx;
                double v_i_sqr = v_x_i*v_x_i;
                // will particle reach boundary?
                future_boundary_bool = (diff_PE > -v_i_sqr);
                if (future_boundary_bool) {
                    xi_f = double(xi_boundary);
                    v_x_f = v_sign * std::sqrt(diff_PE + v_i_sqr);
                    // select method to minimize error in del_tau
                    if (std::abs(v_x_f - v_x_i) > std::abs(v_x_f + v_x_i)) {
                        del_tau = (v_x_f - v_x_i)/accel;
                    } else {
                        del_tau = 2.0 * (xi_f - xi_i)/ (v_x_f + v_x_i) * dx;
                    }
                } else if (!in_cell_bool) {
                    future_boundary_bool = true;
                    xi_boundary = 2 * xi_cell + 1 - xi_boundary;
                    xi_f = double(xi_boundary);
                    diff_PE = 2.0 * accel * (xi_f - xi_i) * dx;
                    v_x_f = - v_sign * std::sqrt(diff_PE + v_i_sqr);
                    del_tau = (v_x_f - v_x_i)/accel;
                }

            }
            if (future_boundary_bool) {
                v_sign = (v_x_f > 0) - (v_x_f < 0);
                xi_cell = xi_cell + v_sign;
                if (xi_boundary == 0) {
                    switch (left_boundary){
                        case 1:
                        case 4:
                            this->energy_loss[thread_id][0] += v_x_f*v_x_f + v_y*v_y + v_z*v_z;
                            this->wall_loss[thread_id][0]++;
                            this->momentum_loss[thread_id][0][0] += v_x_f;
                            this->momentum_loss[thread_id][0][1] += v_y;
                            this->momentum_loss[thread_id][0][2] += v_z;
                            del_part = true;
                            break;
                        case 2:
                            xi_cell = 0;
                            v_x_f = - v_x_f;
                            v_sign = -v_sign;
                            break;
                        case 3:
                            xi_f = double(number_cells);
                            xi_cell = number_cells-1;
                            break;
                    }
                } else if (xi_boundary == number_cells) {
                    switch (right_boundary){
                        case 1:
                        case 4:
                            this->energy_loss[thread_id][1] += v_x_f*v_x_f + v_y*v_y + v_z*v_z;
                            this->wall_loss[thread_id][1]++;
                            this->momentum_loss[thread_id][1][0] += v_x_f;
                            this->momentum_loss[thread_id][1][1] += v_y;
                            this->momentum_loss[thread_id][1][2] += v_z;
                            del_part = true;
                            break;
                        case 2:
                            xi_cell = number_cells-1;
                            v_x_f = -v_x_f;
                            v_sign = -v_sign;
                            break;
                        case 3:
                            xi_f = 0.0;
                            xi_cell = 0;
                            break;
                    }
    
                }
                if (del_part) {
                    // particle is lost, no need to continue
                    break;
                }
            }
            
            time_passed += del_tau;
            del_tau = del_t - time_passed;
            xi_i = xi_f;
            v_x_i = v_x_f;
            
        }
        if (!del_part) {
            size_t new_idx = part_indx-space_delete;
            xi_local[new_idx] = xi_f;
            v_x_local[new_idx] = v_x_f;
            if (use_vy) {
                this->v_y[thread_id][new_idx] = v_y;
            }
            if (use_vz){
                this->v_z[thread_id][new_idx] = v_z;
            }
            if (use_y) {
                this->y[thread_id][new_idx] = this->y[thread_id][part_indx];
            }
            if (use_z) {
                this->z[thread_id][new_idx] = this->z[thread_id][part_indx];
            }
        } else {
            space_delete++;
        }
        
    }
    for (int inj_indx = 0; inj_indx < this->number_unique_injections; inj_indx++) {
        size_t start_indx = last_idx;
        const size_t number_inject_local = this->number_particles_injected[thread_id][inj_indx];
        const std::vector<double>& del_t_array = this->time_step_injected[thread_id][inj_indx];
        last_idx = start_indx + number_inject_local;
        double del_t_local;
        for (size_t part_indx = start_indx; part_indx < last_idx; part_indx++){
            xi_i = xi_local[part_indx];
            v_x_i = v_x_local[part_indx];
            del_t_local = del_t_array[part_indx-start_indx];
            v_sign = (v_x_i > 0) - (v_x_i < 0);
            xi_cell = int(xi_i + v_sign*1e-12); // in case lands on exact boundary, has happened before...
            del_tau = del_t_local;
            time_passed = 0.0;
            del_part = false;
            if (use_vy) {
                v_y = this->v_y[thread_id][part_indx];
            }
            if (use_vz){
                v_z = this->v_z[thread_id][part_indx];
            }
            while (del_tau > t_tol) {
                number_sub_steps++;
                E_field_local = E_field[xi_cell];
                dx = dx_dxi[xi_cell];
                accel = q_over_m * E_field_local;
                v_x_f = v_x_i + accel * del_tau;
                xi_f = xi_i + 0.5 * (v_x_i + v_x_f) * del_tau / dx;
                
                // check if particle outside cell or flipped direction
                in_cell_bool = (int(xi_f+1) == xi_cell+1);
                equal_v_sign_bool = ((v_x_i > 0) == (v_x_f > 0));
                future_boundary_bool = (!in_cell_bool) || (!equal_v_sign_bool);

                if (future_boundary_bool) {
                    
                    if (v_x_i > 0) {
                        xi_boundary = xi_cell + 1;
                    } else {
                        xi_boundary = xi_cell;
                    }
                    double diff_PE = 2.0 * accel * (double(xi_boundary) - xi_i) * dx;
                    double v_i_sqr = v_x_i*v_x_i;
                    // will particle reach boundary?
                    future_boundary_bool = (diff_PE > -v_i_sqr);
                    if (future_boundary_bool) {
                        xi_f = double(xi_boundary);
                        v_x_f = v_sign * std::sqrt(diff_PE + v_i_sqr);
                        // select method to minimize error in del_tau
                        if (std::abs(v_x_f - v_x_i) > std::abs(v_x_f + v_x_i)) {
                            del_tau = (v_x_f - v_x_i)/accel;
                        } else {
                            del_tau = 2.0 * (xi_f - xi_i)/ (v_x_f + v_x_i) * dx;
                        }
                    } else if (!in_cell_bool) {
                        future_boundary_bool = true;
                        xi_boundary = 2 * xi_cell + 1 - xi_boundary;
                        xi_f = double(xi_boundary);
                        diff_PE = 2.0 * accel * (xi_f - xi_i) * dx;
                        v_x_f = - v_sign * std::sqrt(diff_PE + v_i_sqr);
                        del_tau = (v_x_f - v_x_i)/accel;
                    }

                }
                if (future_boundary_bool) {
                    v_sign = (v_x_f > 0) - (v_x_f < 0);
                    xi_cell = xi_cell + v_sign;
                    if (xi_boundary == 0) {
                        switch (left_boundary){
                            case 1:
                            case 4:
                                this->energy_loss[thread_id][0] += v_x_f*v_x_f + v_y*v_y + v_z*v_z;
                                this->wall_loss[thread_id][0]++;
                                this->momentum_loss[thread_id][0][0] += v_x_f;
                                this->momentum_loss[thread_id][0][1] += v_y;
                                this->momentum_loss[thread_id][0][2] += v_z;
                                del_part = true;
                                break;
                            case 2:
                                xi_cell = 0;
                                v_x_f = - v_x_f;
                                v_sign = -v_sign;
                                break;
                            case 3:
                                xi_f = double(number_cells);
                                xi_cell = number_cells-1;
                                break;
                        }
                    } else if (xi_boundary == number_cells) {
                        switch (right_boundary){
                            case 1:
                            case 4:
                                this->energy_loss[thread_id][1] += v_x_f*v_x_f + v_y*v_y + v_z*v_z;
                                this->wall_loss[thread_id][1]++;
                                this->momentum_loss[thread_id][1][0] += v_x_f;
                                this->momentum_loss[thread_id][1][1] += v_y;
                                this->momentum_loss[thread_id][1][2] += v_z;
                                del_part = true;
                                break;
                            case 2:
                                xi_cell = number_cells-1;
                                v_x_f = -v_x_f;
                                v_sign = -v_sign;
                                break;
                            case 3:
                                xi_f = 0.0;
                                xi_cell = 0;
                                break;
                        }
        
                    }
                    if (del_part) {
                        // particle is lost, no need to continue
                        break;
                    }
                }
                
                time_passed += del_tau;
                del_tau = del_t_local - time_passed;
                xi_i = xi_f;
                v_x_i = v_x_f;
                
            }
            if (!del_part) {
                size_t new_idx = part_indx-space_delete;
                xi_local[new_idx] = xi_f;
                v_x_local[new_idx] = v_x_f;
                if (use_vy) {
                    this->v_y[thread_id][new_idx] = v_y;
                }
                if (use_vz){
                    this->v_z[thread_id][new_idx] = v_z;
                }
                if (use_y) {
                    this->y[thread_id][new_idx] = this->y[thread_id][part_indx];
                }
                if (use_z) {
                    this->z[thread_id][new_idx] = this->z[thread_id][part_indx];
                }
            } else {
                space_delete++;
            }
        }
    }
    this->number_particles[thread_id][0] = (last_idx - space_delete);
    this->number_collidable_particles[thread_id][0] = (last_idx - space_delete);
}
