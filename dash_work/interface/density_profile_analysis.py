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
density_dat = pd.read_csv('densities_water.txt', sep='\t', error_bad_lines=False, header=None)
density_dat = np.asarray(density_dat)

# Compute the average for each bin
average_profile = np.mean(density_dat, axis=0)

# Specify data
xdata = np.linspace(0, 1, 200)
ydata = average_profile


# Define the error function profile we fit to
def rho_W(z, rho, h, w_c):

	return 0.5*rho - 0.5*rho*scipy.special.erf((z - h)/(sqrt(2)*w_c))

popt, pcov = curve_fit(rho_W, xdata, ydata)

print(popt)

w_c = 110*popt[2]

# Plot the data
plt.plot(xdata, ydata, label='data')
plt.plot(xdata, rho_W(xdata, *popt), 'g--', label='fit: density=%5.3f g/${cm}^3$, interfacial width=%5.3f $\AA$' % (popt[0], w_c))
plt.xlabel('Normalized box length in z-direction ($L_z$ = 110 $\AA$)')
plt.ylabel('Density (g/${cm}^3$)')
plt.legend()
plt.title('Water density profile (averaged over 1 ps)')
plt.savefig('water_density_profile_fit.png')






