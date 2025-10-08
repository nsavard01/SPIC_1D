# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 16:39:21 2023

@author: Nicolas
"""
from import_libraries_constants import *
from dataSet_class import *


# ---------------------------- Averages -------------------------------------
def plot_ave_density(dataSet, name = "", label = "", marker = 'o', linestyle = '--'):
    if name == "":
        colors = ['b', 'r', 'g', 'k', 'c', 'm', 'y']
        for i,name in enumerate(dataSet.particles.keys()):
            n = dataSet.get_ave_density(name)
            grid = dataSet.get_grid()
            plt.plot(grid, n,  linewidth = 2, linestyle = linestyle, marker = marker, markersize = 2,color = colors[i], label = r'$n_{' + name +  '}$')
        plt.xlabel('Distance (m)')
        plt.ylabel('Particle Density (1/m^3)')
        plt.xlim([grid.min(), grid.max()])
        plt.legend(loc = 'best')
    else:
        if name not in dataSet.particles.keys():
            raise Warning("For average density, particle", name, "does not exist in the dataSet!")
        else:
            n = dataSet.get_ave_density(name)
            grid = dataSet.get_grid()
            plt.plot(grid, n,  linestyle = linestyle, marker = marker, markersize = 2, label = label)
            plt.xlabel('Distance (m)')
            plt.ylabel(name + ' Density (1/m^3)')
            plt.xlim([grid.min(), grid.max()])
 
def plot_ave_phi(dataSet, label = '', marker = 'o', linestyle = '--'):

    phi = dataSet.get_ave_phi()
    grid = dataSet.get_grid()
    plt.plot(grid, phi, marker = marker, linestyle = linestyle, markersize = 2, label = label)
    plt.xlabel('Distance (m)')
    plt.ylabel('Potential (V)')
    plt.xlim([grid.min(), grid.max()])


def plot_ave_temp(dataSet, name = "", label = "", marker = 'o', linestyle = '--'):
    if name == "":
        colors = ['b', 'r', 'g', 'k', 'c', 'm', 'y']
        for i,name in enumerate(dataSet.particles.keys()):
            n = dataSet.get_ave_temp(name)
            grid = dataSet.get_half_grid()
            plt.plot(grid, n,  linewidth = 2, linestyle = linestyle, marker = marker, markersize = 2,color = colors[i], label = r'$n_{' + name +  '}$')
        plt.xlabel('Distance (m)')
        plt.ylabel(name + ' temperature (eV)')
        plt.xlim([0, dataSet.domain_length])
        plt.legend(loc = 'best')
    else:
        if name not in dataSet.particles.keys():
            raise Warning("For average density, particle", name, "does not exist in the dataSet!")
        else:
            n = dataSet.get_ave_temp(name)
            grid = dataSet.get_half_grid()
            plt.plot(grid, n,  linestyle = linestyle, marker = marker, markersize = 2, label = label)
            plt.xlabel('Distance (m)')
            plt.ylabel(name + ' temperature (eV)')
            plt.xlim([0, dataSet.domain_length])


def plot_ave_EDF(dataSet, name = "", label = "", marker = 'o', linestyle = '--'):
    if name == "":
        colors = ['b', 'r', 'g', 'k', 'c', 'm', 'y']
        for i,name in enumerate(dataSet.particles.keys()):
            E_bins, E_counts, E_sizes = dataSet.get_ave_EDF(name)
            norm = np.sum(E_counts * E_sizes)
            E_hist = E_counts / norm
            plt.plot(E_bins, E_hist,  linewidth = 2, linestyle = linestyle, marker = marker, markersize = 2,color = colors[i], label = name)
        plt.xlabel('Energy (eV)')
        plt.ylabel(name + r' EDF (eV$^{-1}$)')
        plt.legend(loc = 'best')
    else:
        if name not in dataSet.particles.keys():
            raise Warning("For average density, particle", name, "does not exist in the dataSet!")
        else:
            E_bins, E_counts, E_sizes = dataSet.get_ave_EDF(name)
            norm = np.sum(E_counts * E_sizes)
            E_hist = E_counts / norm
            plt.plot(E_bins, E_hist, linewidth=2, linestyle=linestyle, marker=marker, markersize=2)
            plt.xlabel('Energy (eV)')
            plt.ylabel(name + r' EDF (eV$^{-1}$)')



# ----------------------- time dependent -----------------------------------

def update_plot(i, ax, x, y_set, x_label, y_label, legend_label):

    ax.clear()
    for i in nameList:
        n = dataSet.getDensity(name, i)
        ax.plot(dataSet.grid, n, 'o-', label = name)
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Particle Density (1/m^3)')
    ax.set_xlim([dataSet.x_min, dataSet.x_max])
    plt.legend(loc = 'lower center')

def animation(dataSet, type, nameList = [],boolMakeAnimation = False, savePath = "temp_fig.gif", pauseTime = 0.05):
    if not nameList:
        nameList = list(dataSet.particles.keys())
    x_grid = dataSet.get_grid()
    x_half_grid = dataSet.get_half_grid()
    x_max = x_grid[-1]
    x_min = x_grid[0]
    # if boolMakeAnimation:
    #     numframes = dataSet.numDiag
    #     fig, ax = plt.subplots()
    #     for name in nameList:
    #         n = dataSet.getDensity(name, 0)
    #         ax.plot(dataSet.grid, n, 'o-', label = name)
    #
    #     ax.set_xlabel('Distance (m)')
    #     ax.set_ylabel(r'Density (m$^{-3}$)')
    #     ax.set_xlim([dataSet.x_min, dataSet.x_max])
    #     #ax.set_ylim([-20, 40])
    #     ani = animation.FuncAnimation(fig, update_plot_Density, frames=range(numframes), interval = 100,fargs=(dataSet, ax, nameList))
    #
    #     ani.save(savePath)
    #     plt.show()



    plt.figure(figsize = (5,4), dpi = 80)
    for u in range(dataSet.num_diag):

        plt.cla()
        if (type == 'phi'):
            y = dataSet.get_phi(u)
            plt.plot(x_grid, y, 'o-')
            plt.xlabel('Distance (m)')
            plt.ylabel(r'Voltage (V)')
        elif (type == 'density'):
            for name in nameList:
                y = dataSet.get_density(name,u)
                plt.plot(x_grid, y, 'o-', label = name)
                plt.xlabel('Distance (m)')
                plt.ylabel(r'Density (1/m^3)')
                plt.legend(loc = 'lower center')
        else:
            for name in nameList:
                y = dataSet.get_temp(name,u)
                plt.plot(x_half_grid, y, 'o-', label = name)
                plt.xlabel('Distance (m)')
                plt.ylabel(r'Temperature (eV)')
                plt.legend(loc = 'lower center')

        plt.xlim([x_min, x_max])
        plt.pause(pauseTime)

# def maxwellEDVF(x, T):
#     return np.sqrt(m_e/2/np.pi / e/ T) * np.exp(- m_e * x**2 / 2 / e/ T)
#
# def plotAveVDF(dataSet, name, CurveFit = False, marker = '', linestyle = '--', label = ''):
#     VHist, Vedge = dataSet.getAveVDF(name)
#     dv = Vedge[1] - Vedge[0]
#     Norm = np.sum(VHist[1:-1] * dv) + 0.5 * (VHist[0] + VHist[-1]) * dv
#     VHist = VHist/Norm
#     plt.plot(Vedge, VHist, marker = marker, linestyle = linestyle, linewidth = 2, label=label)
#     plt.ylabel('VDF (s/m)')
#     plt.xlabel('Velocity (m/s)')
#     if CurveFit:
#         popt, pcov = opt.curve_fit(maxwellEDVF, Vedge, VHist, p0 = [dataSet.T_e])
#         plt.plot(Vedge, maxwellEDVF(Vedge, popt[0]), label = label + 'Fit T_e = ' + '{:2.2f}'.format(popt[0]))
#         plt.legend(loc = 'best')
#
# def plotAveEDF(dataSet, name, CurveFit = False, label = ''):
#     EHist, Ebin = dataSet.getAveEDF(name)
#     EHist = EHist/Ebin
#     Norm = np.trapz(EHist, Ebin)
#     EHist = EHist/Norm
#
#     if CurveFit:
#         f_log = np.log(EHist/np.sqrt(Ebin))
#         fit = np.polyfit(Ebin, f_log, 1)
#         plt.plot(Ebin, EHist/np.sqrt(Ebin), '--', label=label)
#         dist = np.poly1d(fit)
#         plt.plot(Ebin, np.exp(dist(Ebin)), label = label + 'Fit T_e = ' + '{:2.2f}'.format(-1/fit[0]))
#         plt.ylabel(r'EPDF eV$^{-3/2}$')
#         plt.yscale('log')
#         plt.xlabel('Particle Energy (eV)')
#         plt.legend(loc = 'best')
#     else:
#         plt.plot(Ebin, EHist, 'o-', label=label)
#         plt.ylabel('EDF (1/eV)')
#         plt.xlabel('Particle Energy (eV)')
#
#
# def plotWallCurrent(data_point_array, particle_name, dataset_list, label = ""):
#     if (data_point_array.size != len(dataset_list)):
#         raise Exception("Error in size data array vs number of data sets")
#
#
#     val_array = np.zeros(data_point_array.size)
#     for j in range(data_point_array.size):
#         dataset = dataset_list[j]
#         val_array[j] = (dataset.particles[particle_name]['aveDiag']['leftCurrLoss'].values[0] + dataset.particles[particle_name]['aveDiag']['rightCurrLoss'].values[0])
#
#     plt.plot(data_point_array, val_array, 'o-', label = label)
#
#
# def plotAveEPF(dataSet, name, lower_limit = 0.2, label = '', marker = 'o', linestyle = '--'):
#     EHist, Ebin = dataSet.getAveEDF(name)
#     EHist = EHist/Ebin
#     Norm = np.trapz(EHist, Ebin)
#     EHist = EHist/Norm
#     lower_indx = np.where(Ebin < lower_limit)[0][-1]
#     plt.plot(Ebin[lower_indx::], EHist[lower_indx::]/np.sqrt(Ebin[lower_indx::]), linestyle = linestyle, marker = marker, markersize = 2, label=label)
#     plt.ylabel(r'EPDF eV$^{-3/2}$')
#     plt.yscale('log')
#     plt.xlabel('Particle Energy (eV)')
#
# def getAveDensityFiles(dataName, name = 'e', diagNumber = 0):
#     # Given name of file which will have _1, _2, etc appended
#     if ('Explicit' in dataName):
#         data = dataSetExplicit(dataName + '/')
#     else:
#         data = dataSet(dataName + '/')
#     density = data.getDensity(name, diagNumber)
#     fileList = glob.glob(dataName + '_[0-9]')
#     number = len(fileList)
#     for filename in fileList:
#         if ('Explicit' in dataName):
#             data = dataSetExplicit(filename + '/')
#         else:
#             data = dataSet(filename + '/')
#         density = density + data.getDensity(name, diagNumber)
#     print('Finished averaging for', dataName)
#     print('Contained', number+1, 'files')
#     return density/(number+1)
#
# def getAveTempFiles(dataName, name = 'e', diagNumber = 0):
#     # Given name of file which will have _1, _2, etc appended
#     if ('Explicit' in dataName):
#         data = dataSetExplicit(dataName + '/')
#     else:
#         data = dataSet(dataName + '/')
#     temp = data.getTemp(name, diagNumber)
#     fileList = glob.glob(dataName + '_[0-9]')
#     number = len(fileList)
#     for filename in fileList:
#         if ('Explicit' in dataName):
#             data = dataSetExplicit(filename + '/')
#         else:
#             data = dataSet(filename + '/')
#         temp = temp + data.getTemp(name, diagNumber)
#     print('Finished averaging for', dataName)
#     print('Contained', number+1, 'files')
#     return temp/(number+1)
#
# def getAveTempFiles_phaseSpace(dataName, name = 'e'):
#     # Given name of file which will have _1, _2, etc appended
#     if ('Explicit' in dataName):
#         data = dataSetExplicit(dataName + '/')
#     else:
#         data = dataSet(dataName + '/')
#     if (data.scheme == 'CIC'):
#         n_cell = data.Nx
#     else:
#         n_cell = data.Nx-1
#     phase_space = data.getPhaseSpace(name)
#
#     v2_binned_total = np.zeros(n_cell)
#     count_total = np.zeros(n_cell)
#
#     # Assume data is your (N_p, 4) array
#     positions = phase_space[:, 0] - 1 # convert to start with 0
#     velocities = phase_space[:, 1:4]
#
#     # Step 1: Get cell indices (e.g., by flooring to nearest lower int)
#     cell_indices = np.floor(positions).astype(int)
#
#     # Step 2: Compute v^2 = vx^2 + vy^2 + vz^2
#     v_squared = np.sum(velocities ** 2, axis=1)
#
#     # Step 3: Bin v^2 into cell indices using np.bincount
#     # Find max cell index to pre-allocate array of correct size
#
#     # Bin the v^2 values
#     v2_binned_total = v2_binned_total +  np.bincount(cell_indices, weights=v_squared, minlength=n_cell)
#
#     count_total = count_total + np.bincount(cell_indices, minlength=n_cell)
#     fileList = glob.glob(dataName + '_[0-9]')
#     number = len(fileList)
#     for filename in fileList:
#         if ('Explicit' in dataName):
#             data = dataSetExplicit(filename + '/')
#         else:
#             data = dataSet(filename + '/')
#         phase_space = data.getPhaseSpace(name)
#         positions = phase_space[:, 0] - 1  # convert to start with 0
#         velocities = phase_space[:, 1:4]
#
#         # Step 1: Get cell indices (e.g., by flooring to nearest lower int)
#         cell_indices = np.floor(positions).astype(int)
#
#         # Step 2: Compute v^2 = vx^2 + vy^2 + vz^2
#         v_squared = np.sum(velocities ** 2, axis=1)
#
#         # Step 3: Bin v^2 into cell indices using np.bincount
#         # Find max cell index to pre-allocate array of correct size
#
#         # Bin the v^2 values
#         v2_binned_total = v2_binned_total + np.bincount(cell_indices, weights=v_squared, minlength=n_cell)
#
#         count_total = count_total + np.bincount(cell_indices, minlength=n_cell)
#
#     temp = data.particles[name]['mass'] * v2_binned_total / count_total / 3.0 / e
#     print('Finished averaging for', dataName)
#     print('Contained', number+1, 'files')
#     return temp
#
#
# def getAvePhiFiles(dataName, diagNumber = 0):
#     # Given name of file which will have _1, _2, etc appended
#     if ('Explicit' in dataName):
#         data = dataSetExplicit(dataName + '/')
#     else:
#         data = dataSet(dataName + '/')
#     density = data.getPhi(diagNumber)
#     fileList = glob.glob(dataName + '_[0-9]')
#     number = len(fileList)
#     for filename in fileList:
#         if ('Explicit' in dataName):
#             data = dataSetExplicit(filename + '/')
#         else:
#             data = dataSet(filename + '/')
#         density = density + data.getPhi(diagNumber)
#     print('Finished averaging for', dataName)
#     print('Contained', number+1, 'files')
#     return density/(number+1)
#
# def getAveEFieldFiles(dataName, diagNumber = 0):
#     # Given name of file which will have _1, _2, etc appended
#     if ('Explicit' in dataName):
#         data = dataSetExplicit(dataName + '/')
#     else:
#         data = dataSet(dataName + '/')
#     density = data.getEField(diagNumber)[1]
#     grid = data.getEField(diagNumber)[0]
#     fileList = glob.glob(dataName + '_[0-9]')
#     number = len(fileList)
#     for filename in fileList:
#         if ('Explicit' in dataName):
#             data = dataSetExplicit(filename + '/')
#         else:
#             data = dataSet(filename + '/')
#         density = density + data.getEField(diagNumber)[1]
#     print('Finished averaging for', dataName)
#     print('Contained', number+1, 'files')
#     return (grid, density/(number+1))
#
# def get_rms_percentage_ave_density(data, refData, name):
#     ref_grid = refData.grid
#     ref_n = refData.getAveDensity(name)
#
#     grid = data.grid
#     n = data.getAveDensity(name)
#     res = 0.0
#     for u in range(data.Nx):
#         res += ((n[u] - np.interp(grid[u], ref_grid, ref_n))/np.interp(grid[u], ref_grid, ref_n))**2
#     res = np.sqrt(res/data.Nx)
#     print('res is', res)
#
# def get_rms_percentage_ave_phi(data, refData):
#     ref_grid = refData.grid
#     ref_n = refData.getAvePhi()
#
#     grid = data.grid
#     n = data.getAvePhi()
#     res = 0.0
#     for u in range(1,data.Nx-1):
#         res += ((n[u] - np.interp(grid[u], ref_grid, ref_n))/np.interp(grid[u], ref_grid, ref_n))**2
#     res = np.sqrt(res/(data.Nx-2))
#     print('res is', res)
#
#
# def plotAveDensityVsRef(refData, dataList, name, labelList = [''], marker = 'o', linestyle = '--'):
#     #Plotting average density vs a reference density
#     refDensity = refData.getAveDensity(name)
#     refGrid = refData.grid
#     for indx,data in enumerate(dataList):
#         density = data.getAveDensity(name)
#         interpDensity = np.interp(data.grid, refGrid, refDensity)
#         relDensity = (density - interpDensity)
#         if (labelList[0] == ''):
#             labelStr = ''
#         else:
#             labelStr = labelList[indx]
#         plt.plot(data.grid, relDensity, linestyle=linestyle, marker=marker, markersize=2,
#              label= labelStr)
#
# def plotAveDensityVsRef_percent(refData, dataList, name, labelList = [''], marker = 'o', linestyle = '--'):
#     #Plotting average density vs a reference density
#     refDensity = refData.getAveDensity(name)
#     refGrid = refData.grid
#     for indx,data in enumerate(dataList):
#         density = data.getAveDensity(name)
#         interpDensity = np.interp(data.grid, refGrid, refDensity)
#         relDensity = (density - interpDensity)/interpDensity
#         if (labelList[0] == ''):
#             labelStr = ''
#         else:
#             labelStr = labelList[indx]
#         plt.plot(data.grid, relDensity, linestyle=linestyle, marker=marker, markersize=2,
#              label= labelStr)
#
# #-------------------------- Time Dependent ---------------------
#
# def update_plot_Phi(i, dataSet, ax, ylim):
#
#     ax.clear()
#     phi = dataSet.getPhi(i)
#     ax.plot(dataSet.grid, phi, 'o-')
#
#     ax.set_xlabel('Distance (m)')
#     ax.set_ylabel('Potential (V)')
#     ax.set_xlim([dataSet.x_min, dataSet.x_max])
#     if (ylim[0] != 0 and ylim[1] != 0):
#         ax.set_ylim(ylim)
#
#
# def phiAnimation(dataSet, boolMakeAnimation = False, savePath = "Figures/BoundPlasmaPhi.gif", pauseTime = 0.05, ylim = [0,0]):
#
#     if boolMakeAnimation:
#         numframes = dataSet.numDiag
#         fig, ax = plt.subplots()
#         phi = dataSet.getPhi(0)
#         ax.plot(dataSet.grid, phi, 'o-')
#
#         ax.set_xlabel('Distance (m)')
#         ax.set_ylabel('Potential (V)')
#         ax.set_xlim([dataSet.x_min, dataSet.x_max])
#         if (ylim[0] != 0 and ylim[1] != 0):
#             ax.set_ylim(ylim)
#         ani = animation.FuncAnimation(fig, update_plot_Phi, frames=range(numframes), interval = 100,fargs=(dataSet, ax, ylim))
#
#         ani.save(savePath)
#         plt.show()
#
#
#     else:
#         plt.figure(figsize = (5,4), dpi = 80)
#         for y in range(dataSet.numDiag):
#
#             plt.cla()
#
#             phi = dataSet.getPhi(y)
#             plt.plot(dataSet.grid, phi, 'o-')
#             plt.xlabel('Distance (m)')
#             plt.ylabel('Potential (V)')
#             plt.xlim([dataSet.x_min, dataSet.x_max])
#             if (ylim[0] != 0 and ylim[1] != 0):
#                 plt.ylim(ylim)
#             plt.pause(pauseTime)
#

#
# def aveQuantity_vs_parameter(quantity, list_data, parameter_list, label_list, name = 'e'):
#     size_param = len(parameter_list)
#     size_list = len(list_data)
#     for j in range(size_list):
#         ave = np.zeros(size_param)
#         for i in range(size_param):
#             dataSet = list_data[j][i]
#             if (quantity == 'phi'):
#                 q = dataSet.getAvePhi()
#                 q_temp = np.trapz(q, dataSet.grid)
#                 q_temp = q_temp / (dataSet.grid[-1] - dataSet.grid[0])
#             elif (quantity == 'density'):
#                 q = dataSet.getAveDensity(name)
#                 q_temp = np.trapz(q, dataSet.grid)
#                 q_temp = q_temp / (dataSet.grid[-1] - dataSet.grid[0])
#             elif (quantity == 'current'):
#                 q_temp = dataSet.particles[name]['aveDiag']['leftCurrLoss'].values[0] + dataSet.particles[name]['aveDiag']['rightCurrLoss'].values[0]
#             elif (quantity == 'temperature'):
#                 q_temp = dataSet.getAveTemp(name)
#
#             ave[i] = q_temp
#         plt.plot(parameter_list, ave, '--o', label = label_list[j])
#
# def density_std_plots(dataSet, name, amount, label = ''):
#     num_nodes = dataSet.Nx
#     n = np.zeros((amount, num_nodes))
#     num_diag = 0
#     for y in range(dataSet.numDiag-1, dataSet.numDiag - amount-1, -1):
#         n[num_diag,:] = dataSet.getDensity(name, y)
#         num_diag += 1
#
#     # n_tot = np.mean(n, axis = 0)
#     # n_err = np.max(n, axis = 0) - np.min(n,axis = 0)
#     n_err = np.std(n, axis = 0)
#     # plt.errorbar(dataSet.grid, n_tot, yerr = n_err, fmt = '.')
#     plt.plot(dataSet.grid * 1e2, n_err, linewidth = 2, linestyle = '--', marker = 'o', markersize = 2, label = label)
#
# def Efield_std_plots(dataSet, amount, label = ''):
#     E = np.zeros((amount, dataSet.Nx-1))
#     num_diag = 0
#     for y in range(dataSet.numDiag-1, dataSet.numDiag - amount-1, -1):
#         phi = dataSet.getPhi(y)
#         E[num_diag, :] = -np.diff(phi)/np.diff(dataSet.grid)
#         num_diag += 1
#     E_err = np.std(E, axis = 0)
#     plt.plot(0.5 * (dataSet.grid[0:-1] + dataSet.grid[1::]) * 1e2, E_err, linewidth = 2, linestyle = '--', marker = 'o', markersize = 2, label = label)
#
# def temp_plots(dataSet, name, amount):
#     T = np.zeros((amount, dataSet.Nx-1))
#     num_diag = 0
#     for y in range(dataSet.numDiag-1, dataSet.numDiag - amount-1, -1):
#         Temp = dataSet.getTemp(name,y)
#         T[num_diag, :] = Temp
#         num_diag += 1
#     T_err = np.std(T, axis = 0)
#     plt.plot(0.5 * (dataSet.grid[0:-1] + dataSet.grid[1::]), T_err, 'o-', label=y)
#     plt.xlabel('Distance (m)')
#     plt.ylabel(r'Temperature (eV)')
#     plt.xlim([dataSet.x_min, dataSet.x_max])
#
# def plotPhaseSpace(dataSet, name):
#     phaseSpace = dataSet.getPhaseSpace(name)
#     plt.scatter(phaseSpace[:,0], phaseSpace[:,1])
# # def update_plot_PhaseSpace(i, scat, dataSet, name):
# #     phaseSpace = dataSet.getPhaseSpace(name, i)
# #     scat.set_offsets(phaseSpace[:,0:2])
# #     return scat,
# #
# # def PhaseSpaceAnimation(dataSet, name, boolMakeAnimation = False, savePath = 'PhaseSpaceAnimation.gif'):
# #     if name not in dataSet.particles.keys():
# #         raise Warning("Particle", name, "does not exist in the dataSet!")
# #     if boolMakeAnimation:
# #         numframes = dataSet.numDiag
# #         phaseSpace = dataSet.getPhaseSpace(name, 0)
# #         fig = plt.figure()
# #         scat = plt.scatter(phaseSpace[:,0], phaseSpace[:,1])
# #         plt.xlabel('Distance (m)')
# #         plt.ylabel('Speed (m/s)')
# #         vmax = abs(np.max(phaseSpace[:,1]))*2.0
# #         plt.axis([dataSet.grid[0], dataSet.grid[-1], -vmax, vmax])
# #         ani = animation.FuncAnimation(fig, update_plot_PhaseSpace, frames=range(numframes), interval = 100,
# #                                       fargs=(scat, dataSet, name))
# #         ani.save(savePath)
# #         plt.show()
# #     else:
# #         plt.figure(figsize = (5,4), dpi = 80)
# #         for y in range(dataSet.numDiag):
# #             plt.cla()
# #             phaseSpace = dataSet.getPhaseSpace(name, y)
# #             plt.scatter(phaseSpace[:,0], phaseSpace[:,1])
# #             plt.xlabel('Distance (m)')
# #             plt.ylabel('Particle velocity (m/s)')
# #             plt.xlim([0, dataSet.grid[-1]])
# #             plt.pause(0.1)
#
#
#
# #------------ Calculations ---------------------------
#
# def get_thermalization_parameters(dataSet, reaction):
#     nu = 0
#     del_t = dataSet.delT
#     for obj in dataSet.binaryColl[reaction]:
#         nu += obj['aveDiag']['ratio']
#         if (obj['type'] == 'Ionization'):
#             nu_ion = obj['aveDiag']['ratio']/del_t
#     nu = nu / del_t
#     n_e = dataSet.getAveDensity('e').mean()
#
#     EHist, Ebin = dataSet.getAveEDF('e')
#     EHist = EHist/Ebin
#     Norm = np.trapz(EHist, Ebin)
#     EHist = EHist/Norm
#     T_e = np.trapz(EHist*Ebin, Ebin) * 2/3
#     deb = debye_length(T_e, n_e)
#     L = dataSet.grid[-1] - dataSet.grid[0]
#     N_p = dataSet.particles['e']['diag']['N_p'].values[-1]
#     N_D = N_p/(L/deb)
#     print('N_D is', N_D)
#     plas_freq = plasmaFreq(n_e)
#     denom = (N_D**(-2)) + 28 * (1/N_D) * (nu/plas_freq)
#     tau_R = 34.4 / denom / plas_freq
#
#     nu_coulomb = crossSectionCoulomb(T_e, n_e)
#     nu_coulomb = nu_coulomb * n_e * np.sqrt(2 * e * T_e / np.pi / m_e)
#     return tau_R, 1/nu, 1/nu_coulomb, 1/nu_ion
#

