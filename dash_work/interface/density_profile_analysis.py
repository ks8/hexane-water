""" File for plotting average density profiles of interfacial systems, by Kirk Swanson"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys, argparse
from math import *
import numpy as np
import time
import copy, os
import json
import random
import csv
import pandas as pd 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit


# Upload the density profile data
density_dat_water = pd.read_csv('densities_water.txt', sep='\t', error_bad_lines=False, header=None)
density_dat_water = np.asarray(density_dat_water)
density_dat_hexane = pd.read_csv('densities_hexane.txt', sep='\t', error_bad_lines=False, header=None)
density_dat_hexane = np.asarray(density_dat_hexane)

# Compute the average for each bin
average_profile_water = np.mean(density_dat_water, axis=0)
average_profile_hexane = np.mean(density_dat_hexane, axis=0)

# Specify data
xdata = np.linspace(0, 1, 200)
ydata_water = average_profile_water
ydata_hexane = average_profile_hexane


# Define the error function profile we fit to for water
def rho_W(z, rho, h, w_c):

	return 0.5*rho - 0.5*rho*scipy.special.erf((z - h)/(sqrt(2)*w_c))

# Define the error function profile we fit to for hexane
def rho_H(z, rho, h, w_c):

	return 0.5*rho + 0.5*rho*scipy.special.erf((z - h)/(sqrt(2)*w_c))

popt_water, pcov_water = curve_fit(rho_W, xdata, ydata_water)

print(popt_water)

popt_hexane, pcov_hexane = curve_fit(rho_H, xdata, ydata_hexane)

print(popt_hexane)



w_c_water = 110*popt_water[2]
w_c_hexane = 110*popt_hexane[2]

# Plot the data
plt.plot(xdata, ydata_water, label='water data')
plt.plot(xdata, ydata_hexane, label='hexane data')
# plt.plot(xdata, rho_W(xdata, *popt_water), 'g--', label='fit: density=%5.3f g/${cm}^3$, interfacial width=%5.3f $\AA$' % (popt_water[0], w_c_water))
# plt.plot(xdata, rho_H(xdata, *popt_hexane), 'r--', label='fit: density=%5.3f g/${cm}^3$, interfacial width=%5.3f $\AA$' % (popt_hexane[0], w_c_hexane))
plt.plot(xdata, rho_W(xdata, *popt_water), 'g--', label='fit: density=%5.3f g/${cm}^3$' % (popt_water[0]))
plt.plot(xdata, rho_H(xdata, *popt_hexane), 'r--', label='fit: density=%5.3f g/${cm}^3$' % (popt_hexane[0]))
plt.xlabel('Normalized box length in z-direction ($L_z$ = 110 $\AA$)')
plt.ylabel('Density (g/${cm}^3$)')
plt.legend()
plt.title('Water density profile (averaged over 1 ps)')
plt.savefig('interface_density_profile_fit.png')






