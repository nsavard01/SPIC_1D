# -*- coding: utf-8 -*-

"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

import sys
from pathlib import Path
sys.path.append("../SPIC_1D/scripts")
from dataSet_class import *
from plot_production import *
print('Loaded modules')



# MC
data_set = dataSet('./RF_Turner_test/')

# plots
plt.figure()
plot_ave_density(data_set, 'e', label = 'e')
plot_ave_density(data_set, 'He+', label = 'He+')
plt.ylabel(r'Particle density (m$^{-3}$)', fontsize = 14)
plt.xlabel('Distance (m)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower center', fontsize = 12)
plt.tight_layout()
plt.savefig('example_density.pdf')
plt.show()

plt.figure()
plot_ave_phi(data_set)
plt.ylabel(r'Electric Potential (V)', fontsize = 14)
plt.xlabel('Distance (m)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig('example_phi.pdf')
plt.show()

plt.figure()
plot_ave_temp(data_set, 'e', label = 'e')
plot_ave_temp(data_set, 'He+', label = 'He+')
plt.ylabel(r'Particle Temperature (eV)', fontsize = 14)
plt.xlabel('Distance (m)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'upper right', fontsize = 12)
plt.tight_layout()
plt.savefig('example_temp.pdf')
plt.show()
