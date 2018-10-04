"""Script for analysis of q-TIP4P/F simulation output data, by Kirk Swanson"""
# Import modules
import sys
import argparse
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Functions to parse arguments
def create_parser():

	# Create parser and add arguments
	parser = argparse.ArgumentParser(description='Read data files')
	parser.add_argument('-DASH_file', dest='DASH_file', default=None, help='CSV file containing DASH data')
	parser.add_argument('-OpenMM_file', dest='OpenMM_file', default=None, help='CSV file containing OpenMM data')
	parser.add_argument('-filename', dest='filename', default='test', help='Name for output files')

	return parser

# Function to convert arguments
def convert_args(args):

	# Files dictionary
	files={}
	files['DASH_file'] = args.DASH_file
	files['OpenMM_file'] = args.OpenMM_file

	# Options dictionary
	options = {}
	options['filename'] = args.filename

	return files, options

# Function to process the data file
def process_datafile(files, options):

	# Load files and options
	DASH_file = str(files['DASH_file'])
	OpenMM_file = str(files['OpenMM_file'])
	filename = str(options['filename'])

	# Load data file
	DASH_data = np.loadtxt(DASH_file, delimiter=',')
	OpenMM_data = np.loadtxt(OpenMM_file, delimiter=',')

	# Create time sequence
	t = np.arange(0.0, 1000.0, 0.5)

	# Plotting
	plt.figure(figsize=(10, 7.5))
	plt.subplot(2, 1, 1)
	plt.plot(t, DASH_data[:, 2], 'b', label='DASH', zorder=2)
	plt.hlines(np.mean(DASH_data[200:, 2]), 0, 1000, colors='k', linestyles='solid', zorder=3, label='DASH avg', linewidth=2.0)
	plt.plot(t, OpenMM_data[:, 2], 'g', label='OpenMM', zorder=1)
	plt.hlines(np.mean(OpenMM_data[200:, 2]), 0, 1000, colors='k', linestyles='dotted', zorder=3, label='OpenMM avg', linewidth=2.0)
	plt.subplot(2, 1, 1).spines['top'].set_visible(False)	  
	plt.subplot(2, 1, 1).spines['right'].set_visible(False) 
	plt.title('DASH vs OpenMM:'+'\n'+'nBeads=32, rCut=9 $\AA$, step=0.5 fs', fontsize=20, fontweight='bold')
	plt.ylabel('Density (g/mL)', fontsize=12)
	plt.ylim((0.96, 1.1))
	plt.xlim((0, 1000))
	plt.legend(fontsize=10, loc='upper right', bbox_to_anchor=(1.05, 1.3))

	# plt.subplot(3, 1, 2)
	# plt.plot(t, DASH_data[:, 0], 'b', zorder=2)
	# plt.hlines(np.mean(DASH_data[200:, 0]), 0, 1000, colors='k', linestyles='solid', zorder=3)
	# plt.plot(t, OpenMM_data[:, 0], 'g', zorder=1)
	# plt.hlines(np.mean(OpenMM_data[200:, 0]), 0, 1000, colors='k', linestyles='dotted', zorder=3)
	# plt.subplot(3, 1, 2).spines['top'].set_visible(False)	  
	# plt.subplot(3, 1, 2).spines['right'].set_visible(False) 
	# plt.ylabel('Potential Energy (kJ/mol)', fontsize=12)
	# plt.ylim((-44000, -41000))
	# plt.xlim((0, 1000))

	
	plt.subplot(2, 1, 2)
	plt.plot(t, DASH_data[:, 1], 'b', zorder=2)
	plt.hlines(np.mean(DASH_data[200:, 1]), 0, 1000, colors='k', linestyles='solid', zorder=3)
	plt.plot(t, OpenMM_data[:, 1], 'g', zorder=1)
	plt.hlines(np.mean(OpenMM_data[200:, 1]), 0, 1000, colors='k', linestyles='dotted', zorder=3)
	plt.subplot(2, 1, 2).spines['top'].set_visible(False)	  
	plt.subplot(2, 1, 2).spines['right'].set_visible(False) 
	plt.ylabel('Temperature (K)', fontsize=12)
	plt.ylim((8960, 10240))
	plt.xlim((0, 1000))
	plt.xlabel('Simulation Time (ps)', fontsize=12)	

	plt.savefig(filename)


# Main function
def main():

	# Process arguments
	parser = create_parser()
	args = parser.parse_args()
	files, options = convert_args(args)

	# Analyze data
	process_datafile(files, options)

# Run the program
main()


