# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 16:20:31 2023

@author: Nicolas
"""

from import_libraries_constants import *

import tkinter
from tkinter import filedialog
import re



class dataSet:
    
    def __init__(self, filename = ""):
        if filename == "":
            tkinter.Tk().withdraw()
            path = filedialog.askdirectory()
            self.path = path + '/'
        else:
            self.path = filename
        dateTimeBool = os.path.isfile(self.path + 'date_time.dat')
        if (dateTimeBool):
            dateTime = np.loadtxt(self.path + 'date_time.dat', skiprows=1)
            self.dateTime = str(dateTime[0])[0:4] + '-' + str(dateTime[0])[4:6] + '-' + str(dateTime[0])[6:8] + ' ' + \
                            str(int(str(dateTime[1])[0:2]) - 7) + ':' + str(dateTime[1])[2:4] + ':' + str(dateTime[1])[4:6]

        temp_path = self.path + 'initial_condition.dat'
        if (os.path.isfile(temp_path)):
            data = pd.read_csv(temp_path, skiprows = 1, nrows=1, sep='\s+', header = None).values[0,:]
            self.num_mpi = int(data[0])
            self.num_threads = int(data[1])
            scheme = int(data[2])
            if (scheme == 0):
                self.scheme = 'MC-PIC'
            elif (scheme == 1):
                self.scheme = 'EC-PIC'
            elif (scheme == 2):
                self.scheme = 'I-NGP'
            elif (scheme == 3):
                self.scheme = 'I-CIC'
            self.final_sim_time = data[3]
            self.del_t = data[4]
            self.num_diag = int(data[5])
            self.ave_time = data[6]
        temp_path = self.path + '/domain/parameters.dat'
        if (os.path.isfile(temp_path)):
            data = pd.read_csv(temp_path, skiprows=1, nrows=1, sep='\s+', header=None).values[0, :]
            type = int(data[0])
            if (type == 0):
                self.domain_type = 'uniform'
            else:
                self.domain_type = 'non-uniform'
            self.num_nodes = int(data[1])
            self.num_cells = int(data[2])
            h = data[3]
            if (h == 1):
                self.left_boundary = 'Dirichlet'
            elif (h==2):
                self.left_boundary = 'Neumann-Symmetric'
            elif (h == 3):
                self.left_boundary = 'Periodic'
                self.right_boundary = 'Periodic'
            elif (h == 4):
                self.left_boundary = 'RF-Dirichlet'

            h = int(data[4])
            if (h == 1):
                self.right_boundary = 'Dirichlet'
            elif (h==2):
                self.right_boundary = 'Neumann-Symmetric'
            elif (h == 4):
                self.right_boundary = 'RF-Dirichlet'
            self.domain_length = data[5]
        temp_path = self.path + '/phi/parameters.dat'
        if (os.path.isfile(temp_path)):
            data = pd.read_csv(temp_path, skiprows=1, nrows=1, sep='\s+', header=None).values[0, :]
            self.RF_frequency = data[0]/2/np.pi
            self.RF_half_amplitude = data[1]
            self.left_voltage = data[2]
            self.right_voltage = data[3]
            if (self.scheme == 'I-NGP' or self.scheme == 'I-CIC'):
                self.smoothing = (data[4] == 1)
        temp_path = self.path + 'simulation_timing_data.dat'
        if (os.path.isfile(temp_path)):
            data = pd.read_csv(temp_path, skiprows=1, nrows=1, sep='\s+', header=None).values[0, :]
            size = data.size
            if (size == 5):
                self.timing_data = pd.read_csv(temp_path, skiprows = 1, names = ['total', 'field', 'particle', 'null', 'art'], sep='\s+')
        temp_path = self.path + 'global_diagnostic_data.dat'
        if (os.path.isfile(temp_path)):
            data = pd.read_csv(temp_path, skiprows=1, nrows=1, sep='\s+', header=None).values[0, :]
            size = data.size
            if (size == 6):
                self.global_data = pd.read_csv(temp_path, skiprows = 1, names = ['time', 'steps', 'P_x', 'P_y', 'P_z', 'E_tot'], sep='\s+')

        temp_path = self.path + 'field_diagnostics.dat'
        if (os.path.isfile(temp_path)):
            data = pd.read_csv(temp_path, skiprows=1, nrows=1, sep='\s+', header=None).values[0, :]
            size = data.size
            if (size == 2):
                self.field_data = pd.read_csv(temp_path, skiprows=1,
                                               names=['E_tot', 'Gauss Error'], sep='\s+')

        self.particles = {}
        temp_path = self.path + 'charged_particles'
        for filename in os.listdir(temp_path):
            name = filename
            self.particles[name] = {}
            part_path = temp_path + '/' + name + '/'
            part_prop = part_path + 'particle_properties.dat'
            data = pd.read_csv(part_prop, skiprows=1, nrows=1, sep='\s+', header=None).values[0, :]
            self.particles[name]['m'] = data[1]
            self.particles[name]['q'] = data[2]
            self.particles[name]['w'] = data[3]
            self.particles[name]['n_x'] = data[4]
            self.particles[name]['n_v'] = data[5]
            self.particles[name]['diag'] = {}
            part_num = part_path + 'number_diagnostics.dat'
            data = pd.read_csv(part_num, skiprows=1,
                                               names=['N_p', 'left', 'right'], sep='\s+')
            self.particles[name]['diag']['numbers'] = data
            part_P = part_path + 'momentum_diagnostics.dat'
            data = pd.read_csv(part_P, skiprows=1,
                               names=['v_x', 'v_y', 'v_z', 'left_v_x', 'left_v_y', 'left_v_z', 'right_v_x', 'right_v_y', 'right_v_z'], sep='\s+')
            self.particles[name]['diag']['momentum'] = data
            part_E = part_path + 'energy_diagnostics.dat'
            data = pd.read_csv(part_E, skiprows=1,
                               names=['v_sq_x', 'v_sq_y', 'v_sq_z', 'v_sq_tot', 'left_v_sq', 'right_v_sq'], sep='\s+')
            self.particles[name]['diag']['energy'] = data
            if (self.ave_time > 0):
                self.particles[name]['diag']['ave'] = {}
                part_num = part_path + 'number_diagnostics_average.dat'
                data = pd.read_csv(part_num, skiprows=1,
                                   names=['N_p', 'left', 'right'], sep='\s+')
                self.particles[name]['diag']['ave']['numbers'] = data
                part_P = part_path + 'momentum_diagnostics_average.dat'
                data = pd.read_csv(part_P, skiprows=1,
                                   names=['v_x', 'v_y', 'v_z', 'left_v_x', 'left_v_y', 'left_v_z', 'right_v_x',
                                          'right_v_y', 'right_v_z'], sep='\s+')
                self.particles[name]['diag']['ave']['momentum'] = data
                part_E = part_path + 'energy_diagnostics_average.dat'
                data = pd.read_csv(part_E, skiprows=1,
                                   names=['v_sq_x', 'v_sq_y', 'v_sq_z', 'v_sq_tot', 'left_v_sq', 'right_v_sq'],
                                   sep='\s+')
                self.particles[name]['diag']['ave']['energy'] = data
            if (os.path.isdir(part_path + 'null_collision')):
                self.particles[name]['diag']['null'] = {}
                t_path = part_path + 'null_collision'
                for target in os.listdir(t_path):
                    self.particles[name]['diag']['null'][target] = []
                    dir_name = t_path + '/' + target
                    num_coll = int(sum(1 for f in os.listdir(dir_name) if f.startswith("collision_properties_") and os.path.isfile(os.path.join(dir_name, f))))
                    for i in range(num_coll):
                        temp_dict = {}
                        file_prop = dir_name + '/collision_properties_' + str(i) + '.dat'
                        file_diag = dir_name + '/collision_diagnostics_' + str(i) + '.dat'
                        data = pd.read_csv(file_prop, skiprows=1, nrows=1, sep='\s+', header=None).values[0, :]
                        type = int(data[1])
                        if (type == 1):
                            temp_dict['type'] = 'Elastic'
                        elif (type == 2):
                            temp_dict['type'] = 'Ionization'
                        elif (type == 3):
                            temp_dict['type'] = 'Excitation'
                        elif (type == 4):
                            temp_dict['type'] = 'Charge-Exchange'
                        temp_dict['E_thres'] = data[2]
                        temp_dict['max_sig'] = data[3]
                        temp_dict['E_peak'] = data[4]
                        data = pd.read_csv(file_diag,skiprows=1,
                               names=['num_coll', 'num_tot', 'E_loss', 'E_i'], sep='\s+')
                        temp_dict['diag'] = data
                        if (self.ave_time > 0):
                            file_diag = dir_name + '/collision_diagnostics_average_' + str(i) + '.dat'
                            data = pd.read_csv(file_diag, skiprows=1,
                                               names=['num_coll', 'num_tot', 'E_loss', 'E_i'], sep='\s+')
                            temp_dict['ave'] = data
                        self.particles[name]['diag']['null'][target].append(temp_dict)




        self.targets = {}
        temp_path = self.path + 'target_particles'
        for filename in os.listdir(temp_path):
            name = filename
            self.targets[name] = {}
            part_path = temp_path + '/' + name + '/'
            part_prop = part_path + 'particle_properties.dat'
            part_diag = part_path + 'ave_diagnostics.dat'
            data = pd.read_csv(part_prop, skiprows=1, nrows=1, sep='\s+', header=None).values[0, :]
            self.targets[name]['m'] = data[1]
            data = pd.read_csv(part_diag, skiprows=1,
                               names=['n', 'T'], sep='\s+')
            self.targets[name]['diag'] = data

        temp_file = self.path + 'non_linear_solver_properties.dat'
        if (os.path.isfile(temp_file)):
            data = pd.read_csv(temp_file, skiprows=1, nrows=1, sep='\s+', header=None).values[0]
            self.non_linear_type = data[0]
            self.non_linear_eps_a = data[1]
            self.non_linear_eps_r = data[2]
            self.non_linear_max_iter = data[3]
            self.non_linear_params = data[5::]
            self.non_linear_diag = pd.read_csv(self.path + 'non_linear_solver_diagnostics.dat', skiprows=1,
                               names=['time', 'res_norm', 'iterations'], sep='\s+')

        temp_file = self.path + 'global_averaging_diag.dat'
        if (os.path.isfile(temp_file)):
            self.averaging_diag = pd.read_csv(temp_file, skiprows=1,
                                               names=['Ave time', 'sim time', 'EDF time', 'diag time', 'final time', 'res V', 'res n'], sep='\s+')

    def get_grid(self):
        return np.fromfile(self.path + '/domain/grid.dat', dtype = 'float')

    def get_half_grid(self):
        x = np.fromfile(self.path + '/domain/grid.dat', dtype = 'float')
        return 0.5 * (x[0:-1] + x[1::])

    def get_dx(self):
        data = np.fromfile(self.path + '/domain/dx_dxi.dat', dtype = 'float')
        if (self.domain_type == 'uniform'):
            return data[0]
        else:
            return data

    def get_phi(self, diag_number):
        return np.fromfile(self.path + '/phi/potential_' + str(diag_number) +  '.dat', dtype = 'float')

    def get_ave_phi(self):
        if (self.ave_time > 0):
            return np.fromfile(self.path + '/phi/potential_average.dat', dtype = 'float')

    def get_density(self, name, diag_number):
        temp_path = self.path + 'charged_particles/' + name + '/density/density_' + str(diag_number) + '.dat'
        return np.fromfile(temp_path, dtype = 'float')

    def get_ave_density(self, name):
        if (self.ave_time > 0):
            temp_path = self.path + 'charged_particles/' + name + '/density/density_average.dat'
            return np.fromfile(temp_path, dtype = 'float')

    def get_temp(self, name, diag_number):
        temp_path = self.path + 'charged_particles/' + name + '/temperature/cell_temp_' + str(diag_number) + '.dat'
        return np.fromfile(temp_path, dtype = 'float')

    def get_ave_temp(self, name):
        if (self.ave_time > 0):
            temp_path = self.path + 'charged_particles/' + name + '/temperature/cell_temp_average.dat'
            return np.fromfile(temp_path, dtype = 'float')

    def get_ave_EDF(self, name):
        if (self.ave_time > 0):
            temp_path = self.path + 'charged_particles/' + name + '/'
            E_bins = np.fromfile(temp_path + 'EDF_average_bins.dat', dtype = 'float')
            E_counts = np.fromfile(temp_path + 'EDF_average_counts.dat', dtype = 'float')
            E_sizes = np.fromfile(temp_path + 'EDF_average_bin_sizes.dat', dtype = 'float')
            return E_bins, E_counts, E_sizes



    
    
        
            