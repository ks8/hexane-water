"Script to analyze interfacial tension, by Kirk Swanson"
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys, argparse
from math import *
import numpy as np
import time
import copy, os
from ast import literal_eval as make_tuple
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
import math
import csv
import pandas as pd 
import scipy
from scipy.optimize import curve_fit
from scipy import stats

# Function to parse arguments 
def create_parser():

    # Create parser and add arguments
    parser = argparse.ArgumentParser(description='Read data files and options')
    parser.add_argument('-pressure_file', dest='pressure_file', default=None, help='Name of pressure tensor text file')
    parser.add_argument('-z_len', dest='z_len', default=100, help='Name of pressure tensor text file')
    parser.add_argument('-index_start', dest='index_start', default=0, help='Set the initial index for density profile selection')
    parser.add_argument('-index_end', dest='index_end', default=2000, help='Set the final index for density profile selection')

    return parser

# Auxiliary Functions
""" Function to parse args """
def convert_args(args):

    # Files dictionary
    files={}
    files['pressure_file'] = args.pressure_file

  
  	# Options dictionary
    options = {}
    options['z_len'] = args.z_len
    options['index_start'] = args.index_start
    options['index_end'] = args.index_end


    return files, options

# Main processing function
def process_datafile(files, options):

	# Options
	z_len = float(options['z_len'])
	index_start = int(options['index_start'])
	index_end = int(options['index_end'])

	# Upload the density profile data
	tension_data = pd.read_csv(files['pressure_file'], sep='\t', error_bad_lines=False, header=None)
	tension_data = np.asarray(tension_data)

	# Select subset of data to work with
	tension_data = tension_data[index_start:index_end, :]
	num_samples = len(tension_data)



	###################### ANALYSIS #######################################


	average_tension = np.mean(tension_data, axis=0)



	# Multiply Lz to get Angstroms.  Multiply pressure tensor elements to convert from atmospheres to Newtons.  Multiply by 1000 to get milli Newtons. 
	def tension(Pxx, Pyy, Pzz, Lz):

		tension = (Lz*pow(10.0, -10)*0.5*(Pzz*101324.997664 - 0.5*(Pxx*101324.997664 + Pyy*101324.997664)))*1000.0

		return tension

	print(tension(average_tension[0], average_tension[1], average_tension[2], z_len))



	###################### Block Averaging #######################################
	# # Discover the right block size
	# SE_values = []

	# for b in range(1000):

	# 	print(b+1)

	# 	num_SE_blocks = float(num_samples)/float(b+1)
	# 	num_SE_blocks = int(num_SE_blocks)

	# 	# Split the data into the scpecified number of blocks
	# 	tension_blocks = list(np.zeros(num_SE_blocks))

	# 	for i in range(num_SE_blocks):

	# 		tension_blocks[i] = np.mean(tension_data[i*int(num_samples/num_SE_blocks):int((i+1)*num_samples/num_SE_blocks),:], axis=0)

	# 	tension_blocks = np.asarray(tension_blocks)

	# 	# Compute bulk densities and width values for each data block
	# 	interfacial_tensions = list(np.zeros(num_SE_blocks))


	# 	for i in range(num_SE_blocks):



	# 		interfacial_tensions[i] = tension(tension_blocks[i][0], tension_blocks[i][1], tension_blocks[i][2], z_len)

	# 	SE_values.append(list([b, stats.sem(interfacial_tensions)]))

	# print(SE_values)








def main(argv):

    parser = create_parser()
    args = parser.parse_args()
    files, options = convert_args(args)

    process_datafile(files, options)

if __name__ == "__main__":
    main(sys.argv[1:])