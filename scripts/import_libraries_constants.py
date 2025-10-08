# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 16:18:53 2023

@author: Nicolas
"""
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import scipy.constants as constants
import scipy.optimize as opt
import pandas as pd
import scipy
import glob, os
import matplotlib.animation as animation
import math

eps_0 = constants.epsilon_0
c = constants.c
m_e = constants.m_e
m_p = constants.m_p
mu_0 = constants.mu_0
k_boltz = constants.k
e = constants.e
m_amu = 1.66053906660e-27

def debye_length(T_e, n):
    """ Distance of potential drop in plasma, T_e in V"""
    return np.sqrt(eps_0*T_e/n/(e))

def plasmaFreq(n):
    """ Distance of potential drop in plasma, T_e in V"""
    return np.sqrt(n * e**2 / m_e / eps_0)

def gyroradius(voltage, mass, charge, B): #voltage in V
    energy = voltage*charge

    gamma = energy/(mass*(c**2)) + 1
    beta = np.sqrt(1 - 1/(gamma**2))
    
    v = beta*c
    R = gamma*mass*v/(charge*B)
    return R

def crossSectionCoulomb(T_e, n_e):
    b_0 = e / (4 * np.pi * eps_0) / T_e
    omega = 2 * debye_length(T_e, n_e)/b_0
    return (8/np.pi) * b_0**2 * np.log(omega)

def maxwellDistributionE(E, T_e):
    return 2 * np.sqrt(E/np.pi) * (T_e)**(-1.5) * np.exp(-E/T_e)

def maxwellDistributionVx(v_x, T, m):
    v_therm_sqr = 2 * e * T/m
    return 1/np.sqrt(v_therm_sqr * np.pi) * np.exp(-v_x**2 / v_therm_sqr)