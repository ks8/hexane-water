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
from scipy.optimize import leastsq
from scipy import stats

# Function to parse arguments 
def create_parser():

    # Create parser and add arguments
    parser = argparse.ArgumentParser(description='Read data files')
    parser.add_argument('-water_density_dat', dest='water_density_dat', default=None, help='Name of water density profile text file')
    parser.add_argument('-hexane_density_dat', dest='hexane_density_dat', default=None, help='Name of hexane density profile text file')
    parser.add_argument('-nBeads', dest='nBeads', default=1, help='Number of path integral beads')
    parser.add_argument('-z_len', dest='z_len', default=100.0, help='Set the z simulation box length')
    parser.add_argument('-bulk_water_zlo_bin', dest='bulk_water_zlo_bin', default=60, help='Set the lower z bound histogram bin of bulk water')
    parser.add_argument('-bulk_water_zhi_bin', dest='bulk_water_zhi_bin', default=80, help='Set the upper z bound histogram bin of bulk water')
    parser.add_argument('-bulk_hexane_zlo_bin', dest='bulk_hexane_zlo_bin', default=110, help='Set the lower z bound histogram bin of bulk hexane')
    parser.add_argument('-bulk_hexane_zhi_bin', dest='bulk_hexane_zhi_bin', default=124, help='Set the upper z bound histogram bin of bulk hexane')
    parser.add_argument('-num_SE_blocks', dest='num_SE_blocks', default=5, help='Set the number of blocks for SE calculation')
    parser.add_argument('-filename', dest='filename', default='test', help='Set the filename for PNG')
    parser.add_argument('-index_start', dest='index_start', default=0, help='Set the initial index for density profile selection')
    parser.add_argument('-index_end', dest='index_end', default=2000, help='Set the final index for density profile selection')

    return parser

# Auxiliary Functions
""" Function to parse args """
def convert_args(args):

    # Files dictionary
    files={}
    files['water_density_dat'] = args.water_density_dat
    files['hexane_density_dat'] = args.hexane_density_dat
  
  	# Options dictionary
    options = {}
    options['nBeads'] = args.nBeads
    options['z_len'] = args.z_len
    options['bulk_water_zlo_bin'] = args.bulk_water_zlo_bin
    options['bulk_water_zhi_bin'] = args.bulk_water_zhi_bin
    options['bulk_hexane_zlo_bin'] = args.bulk_hexane_zlo_bin
    options['bulk_hexane_zhi_bin'] = args.bulk_hexane_zhi_bin
    options['num_SE_blocks'] = args.num_SE_blocks
    options['filename'] = args.filename
    options['index_start'] = args.index_start
    options['index_end'] = args.index_end

    return files, options

# Main processing function
def process_datafile(files, options):

	# Parse options
	nBeads = int(options['nBeads'])
	box_len_av = float(options['z_len'])
	bulk_water_zlo_bin = int(options['bulk_water_zlo_bin'])
	bulk_water_zhi_bin = int(options['bulk_water_zhi_bin'])
	bulk_hexane_zlo_bin = int(options['bulk_hexane_zlo_bin'])
	bulk_hexane_zhi_bin = int(options['bulk_hexane_zhi_bin'])
	num_SE_blocks = int(options['num_SE_blocks'])
	filename = str(options['filename'])
	index_start = int(options['index_start'])
	index_end = int(options['index_end'])

	# Upload the density profile data
	density_dat_water = pd.read_csv(files['water_density_dat'], sep='\t', error_bad_lines=False, header=None)
	density_dat_water = np.asarray(density_dat_water)
	density_dat_hexane = pd.read_csv(files['hexane_density_dat'], sep='\t', error_bad_lines=False, header=None)
	density_dat_hexane = np.asarray(density_dat_hexane)
	
	# Choose subset of data to work with
	density_dat_water = density_dat_water[index_start:index_end, :]
	density_dat_hexane = density_dat_hexane[index_start:index_end, :]
	num_samples = len(density_dat_water)
	num_bins = len(density_dat_water[0, :])

	# Compute the average for each bin
	average_profile_water = np.mean(density_dat_water, axis=0)
	average_profile_hexane = np.mean(density_dat_hexane, axis=0)

	# Remove double counting of beads 
	average_profile_water = np.multiply(np.array([1/nBeads]), average_profile_water)
	average_profile_hexane = np.multiply(np.array([1/nBeads]), average_profile_hexane)

	# Define global variables
	global rho_w, rho_h

	# Compute the average bulk density for water, using the specified bulk region
	rho_w = np.mean(average_profile_water[bulk_water_zlo_bin:bulk_water_zhi_bin])

	# Compute the average bulk density for hexane, using the specified bulk region
	rho_h = np.mean(average_profile_hexane[bulk_hexane_zlo_bin:bulk_hexane_zhi_bin])

	# Return the bulk densities
	print('bulk water density: '+str(rho_w))
	print('bulk hexane density: '+str(rho_h))

	# Specify data for density profile fitting
	xdata = np.linspace(0, 1, num_bins)
	ydata_water = average_profile_water
	ydata_hexane = average_profile_hexane

	# Define the error function profile we fit to for water
	def rho_W(z, h, w_c):

		return 0.5*rho_w - 0.5*rho_w*scipy.special.erf((z - h)/(sqrt(2)*w_c))

	# Define the error function profile we fit to for hexane
	def rho_H(z, h, w_c):

		return 0.5*rho_h + 0.5*rho_h*scipy.special.erf((z - h)/(sqrt(2)*w_c))

	# Fit the water density profile and hexane density profile, defining the interfacial region as between the specified bulk regions
	popt_water, pcov_water = curve_fit(rho_W, xdata[bulk_water_zhi_bin:bulk_hexane_zlo_bin], ydata_water[bulk_water_zhi_bin:bulk_hexane_zlo_bin])
	popt_hexane, pcov_hexane = curve_fit(rho_H, xdata[bulk_water_zhi_bin:bulk_hexane_zlo_bin], ydata_hexane[bulk_water_zhi_bin:bulk_hexane_zlo_bin])

	# Print the intrinsic and thermal widths
	print('separate fitting intrinsic width: '+str(box_len_av*(popt_hexane[0] - popt_water[0])))
	print('water thermal width: '+str(box_len_av*popt_water[1]))
	print('hexane thermal width: '+str(box_len_av*popt_hexane[1]))

	# Attempt at simultaneous fitting
	def my_residuals(parameters, xdata, ydata_water, ydata_hexane, density_w, density_h):

		# Unpack the parameters from the parameters array
		h_W = parameters[0]
		h_H = parameters[1]
		w_c = parameters[2]

		res1 = ydata_water - (0.5*density_w - 0.5*density_w*scipy.special.erf((xdata - h_W)/(sqrt(2)*w_c)))
		res2 = ydata_hexane - (0.5*density_h + 0.5*density_h*scipy.special.erf((xdata - h_H)/(sqrt(2)*w_c)))

		return np.concatenate((res1, res2))

	par_init = np.array([popt_water[0], popt_hexane[0], popt_water[1]])

	fit_total, cov_total = leastsq(my_residuals, par_init, args=(xdata[bulk_water_zhi_bin:bulk_hexane_zlo_bin], ydata_water[bulk_water_zhi_bin:bulk_hexane_zlo_bin], ydata_hexane[bulk_water_zhi_bin:bulk_hexane_zlo_bin], rho_w, rho_h))

	print('total intrinsic width:  '+str(box_len_av*(fit_total[1] - fit_total[0])))
	print('total thermal width:  '+str(box_len_av*fit_total[2]))

	# SE calculation, split the data into the specified number of blocks 
	means_water = list(np.zeros(num_SE_blocks))
	means_hexane = list(np.zeros(num_SE_blocks))

	for i in range(num_SE_blocks):

		means_water[i] = np.mean(density_dat_water[i*int(num_samples/num_SE_blocks):int((i+1)*num_samples/num_SE_blocks),:], axis=0)
		means_hexane[i] = np.mean(density_dat_hexane[i*int(num_samples/num_SE_blocks):int((i+1)*num_samples/num_SE_blocks),:], axis=0)

		# Normalize for bead number
		means_water[i] = np.multiply(np.array([1/nBeads]), means_water[i])
		means_hexane[i] = np.multiply(np.array([1/nBeads]), means_hexane[i])

	SE_water = list(np.zeros(num_bins))
	SE_hexane = list(np.zeros(num_bins))

	for i in range(num_bins):

		water_list = []
		hexane_list = []

		for j in range(num_SE_blocks):

			water_list.append(means_water[j][i])
			hexane_list.append(means_hexane[j][i])

		SE_water[i] = stats.sem(water_list)
		SE_hexane[i] = stats.sem(hexane_list)

	intrinsic_widths = list(np.zeros(num_SE_blocks))
	thermal_widths = list(np.zeros(num_SE_blocks))
	h_water = list(np.zeros(num_SE_blocks))
	h_hexane = list(np.zeros(num_SE_blocks))
	box_widths = list(np.zeros(num_SE_blocks))

	for i in range(num_SE_blocks):

		fit_test, cov_test = leastsq(my_residuals, par_init, args=(xdata[bulk_water_zhi_bin:bulk_hexane_zlo_bin], means_water[i][bulk_water_zhi_bin:bulk_hexane_zlo_bin], means_hexane[i][bulk_water_zhi_bin:bulk_hexane_zlo_bin], np.mean(means_water[i][bulk_water_zlo_bin:bulk_water_zhi_bin]), np.mean(means_hexane[i][bulk_hexane_zlo_bin:bulk_hexane_zhi_bin])))

		intrinsic_widths[i] = box_len_av*(fit_test[1] - fit_test[0])
		thermal_widths[i] = box_len_av*fit_test[2]

		# For plotting
		h_water[i] = fit_test[0]
		h_hexane[i] = fit_test[1]
		box_widths = fit_test[2]


	print('intrinsic width:  '+str(np.mean(intrinsic_widths))+'  +/- '+str(stats.sem(intrinsic_widths)))
	print('thermal width:  '+str(np.mean(thermal_widths))+'  +/- '+str(stats.sem(thermal_widths)))


	# Plot the data
	plt.plot(xdata, ydata_water, label='water data')
	plt.plot(xdata, ydata_hexane, label='hexane data')

	# # Plot the curve individual fits
	# plt.plot(xdata, rho_W(xdata, *popt_water), 'g--', label='fit: water')
	# plt.plot(xdata, rho_H(xdata, *popt_hexane), 'r--', label='fit: hexane')

	# # Plot the curve fits from simultaneous fit 
	# plt.plot(xdata, rho_W(xdata, fit_total[0], fit_total[2]), 'g--', label='fit: water')
	# plt.plot(xdata, rho_H(xdata, fit_total[1], fit_total[2]), 'r--', label='fit: hexane')

	# Plot the curve fits from the block averaged simultaneous fits
	plt.plot(xdata, rho_W(xdata, np.mean(h_water), np.mean(box_widths)), 'g--', label='fit: water')
	plt.plot(xdata, rho_H(xdata, np.mean(h_hexane), np.mean(box_widths)), 'r--', label='fit: hexane')

	# Plot the standard error measurements
	for i in range(num_bins):
		plt.fill_between(xdata, ydata_water - SE_water, ydata_water + SE_water,color='violet',alpha=0.5)

	for i in range(num_bins):
		plt.fill_between(xdata, ydata_hexane - SE_hexane, ydata_hexane + SE_hexane,color='violet',alpha=0.5)
		
	# Plot labels
	plt.xlabel('Box length in z-direction')
	plt.ylabel('Density (g/${cm}^3$)')
	plt.legend()
	plt.title('Density profile')
	plt.savefig('interface_density_profile_'+filename+'.png')

def main(argv):

    parser = create_parser()
    args = parser.parse_args()
    files, options = convert_args(args)

    process_datafile(files, options)

if __name__ == "__main__":
    main(sys.argv[1:])






