""" File for plotting average orientation profiles of interfacial systems, by Kirk Swanson"""
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
from scipy import stats

# Function to parse arguments for one DASH restart file
def create_parser():

    # Create parser and add arguments
    parser = argparse.ArgumentParser(description='Read data files')
    parser.add_argument('-orientations_dat', dest='orientations_dat', default=None, help='Name of text file containing orientation profile data')
    parser.add_argument('-water_density_dat', dest='water_density_dat', default=None, help='Name of water density profile text file')
    parser.add_argument('-hexane_density_dat', dest='hexane_density_dat', default=None, help='Name of hexane density profile text file')
    parser.add_argument('-zlo_bin', dest='zlo_bin', default=600, help='Set the lower z bound histogram bin')
    parser.add_argument('-zhi_bin', dest='zhi_bin', default=800, help='Set the upper z bound histogram bin')
    parser.add_argument('-index_start', dest='index_start', default=0, help='Set the initial index for density profile selection')
    parser.add_argument('-index_end', dest='index_end', default=2000, help='Set the final index for density profile selection')
    parser.add_argument('-nBeads', dest='nBeads', default=1, help='Number of path integral beads')
    parser.add_argument('-filename', dest='filename', default='test', help='Set the filename for PNG')
    parser.add_argument('-num_SE_blocks', dest='num_SE_blocks', default=5, help='Set the number of blocks for SE calculation')

    return parser

""" Function to parse args """
def convert_args(args):

    # Files dictionary
    files={}
    files['orientations_dat'] = args.orientations_dat
    files['water_density_dat'] = args.water_density_dat
    files['hexane_density_dat'] = args.hexane_density_dat

    options = {}
    options['zlo_bin'] = args.zlo_bin
    options['zhi_bin'] = args.zhi_bin
    options['index_start'] = args.index_start
    options['index_end'] = args.index_end
    options['nBeads'] = args.nBeads
    options['filename'] = args.filename
    options['num_SE_blocks'] = args.num_SE_blocks

    return files, options

# Main processing function
def process_datafile(files, options):

	# Parse options
	zlo_bin = int(options['zlo_bin'])
	zhi_bin = int(options['zhi_bin'])
	index_start = int(options['index_start'])
	index_end = int(options['index_end'])
	nBeads = int(options['nBeads'])
	filename = str(options['filename'])
	num_SE_blocks = int(options['num_SE_blocks'])

	# Upload the orientation profile data
	orientation_dat_water = pd.read_csv(files['orientations_dat'], sep='\t', error_bad_lines=False, header=None)
	orientation_dat_water = np.asarray(orientation_dat_water)
	orientation_dat_water = orientation_dat_water[index_start:index_end, :]
	orientation_dat_water[orientation_dat_water < -1] = np.nan

	# Compute the average for each bin
	average_orientation_water = np.nanmean(orientation_dat_water, axis=0)
	average_orientation_water[0:zlo_bin] = np.nan
	average_orientation_water[zhi_bin:] = np.nan

	# Upload the density profile data
	density_dat_water = pd.read_csv(files['water_density_dat'], sep='\t', error_bad_lines=False, header=None)
	density_dat_water = np.asarray(density_dat_water)
	density_dat_hexane = pd.read_csv(files['hexane_density_dat'], sep='\t', error_bad_lines=False, header=None)
	density_dat_hexane = np.asarray(density_dat_hexane)

	# Choose subset of data to work with 
	density_dat_water = density_dat_water[index_start:index_end, :]
	density_dat_hexane = density_dat_hexane[index_start:index_end, :]

	# Remove double counting of beads
	density_dat_water = np.multiply(np.array(1/nBeads), density_dat_water)
	density_dat_hexane = np.multiply(np.array(1/nBeads), density_dat_hexane)

	# Compute number of bins used in analysis 
	num_samples = len(density_dat_water)
	num_bins = len(density_dat_water[0, :])

	# Compute the average density profile for each bin
	average_profile_water = np.mean(density_dat_water, axis=0)
	average_profile_hexane = np.mean(density_dat_hexane, axis=0)

	# Rescale density profiles for plotting purposes
	average_profile_water = np.multiply(np.array([0.15]), average_profile_water)
	average_profile_water[:] = [x - 0.05 for x in average_profile_water]
	average_profile_hexane = np.multiply(np.array([0.15]), average_profile_hexane)
	average_profile_hexane[:] = [x - 0.05 for x in average_profile_hexane]

	# Split the data into the specified number of blocks and compute average profiles
	orientation_blocks_water = list(np.zeros(num_SE_blocks))

	for i in range(num_SE_blocks):

		orientation_blocks_water[i] = np.nanmean(orientation_dat_water[i*int(num_samples/num_SE_blocks):int((i+1)*num_samples/num_SE_blocks),:], axis=0)
		orientation_blocks_water[i][0:zlo_bin] = np.nan
		orientation_blocks_water[i][zhi_bin:] = np.nan

	# Compute SE measurements for the density profiles for plotting
	SE_orientation = list(np.zeros(num_bins))

	for i in range(num_bins):

		orientation_list = []

		for j in range(num_SE_blocks):

			orientation_list.append(orientation_blocks_water[j][i])

		SE_orientation[i] = stats.sem(orientation_list)


	# Plot the data
	xdata = np.linspace(0, 1, num_bins)
	plt.plot(xdata, average_orientation_water, label='water orientation', color='blue')
	plt.plot(xdata, average_profile_water, label='scaled water density', color='darkgrey')
	plt.plot(xdata, average_profile_hexane, label='scaled hexane density', color='dimgrey')

	# Plot the standard error measurements
	for i in range(num_bins):
		plt.fill_between(xdata, average_orientation_water - SE_orientation, average_orientation_water + SE_orientation,color='violet',alpha=0.5)

	# Plot labels
	plt.xlabel('Normalized box length in z-direction')
	plt.ylabel('$S_{zz}$')
	plt.legend()
	plt.title('Water orientation profile')

	plt.savefig(filename+'.png')

def main(argv):

    parser = create_parser()
    args = parser.parse_args()
    files, options = convert_args(args)

    process_datafile(files, options)

if __name__ == "__main__":
    main(sys.argv[1:])