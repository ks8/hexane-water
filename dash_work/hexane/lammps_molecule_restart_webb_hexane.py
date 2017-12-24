"Script to build a restart LAMMPS file from a DASH restart file of hexane originally read in by the LAMMPS reader, by Kirk Swanson"
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

# Function to parse arguments for LAMMPS and DASH files
def create_parser():

	# Create parser and add arguments
	parser = argparse.ArgumentParser(description='Read data files')
	parser.add_argument('-lammps_file', dest='lammps_file', default=None, help='Name of original LAMMPS input file')
	parser.add_argument('-dash_file', dest='dash_file', default=None, help='Name of DASH restart file')
	parser.add_argument('-out', dest='outfile', default="hexane_restart.txt", help='Desired output file name.  (default: hexane_restart.txt')
	
	return parser

# Function to convert arguments into a dictionary
def convert_args(args):

	# Files dictionary
	files={}
	files['lammps_file'] = args.lammps_file
	files['dash_file'] = args.dash_file
	files['out'] = args.outfile
	
	options = {}

	# Print confirmation
	print("**************************")
	print("# Data file 1: ", files['lammps_file'])
	print("# Data file 2: ", files['dash_file'])
	print("# Final datafile for LAMMPS input ", files['out'])
	print("**************************")
	print(" ")

	return files, options

# Function to find the index location of a substring in a string list
def indices_containing_substring(the_list, substring):
	indices = []
	for i, s in enumerate(the_list):
		if substring in s:
			indices.append(i)	
	return indices

def process_datafile(files):

	# Read the origina lammps file
	f1 = open(files['lammps_file'], 'r')
	contents_lammps = f1.readlines()
	f1.close()

	# Read the DASH restart file
	f2 = open(files['dash_file'], 'r')
	contents_dash = f2.readlines()
	f2.close()

	# Calculate number of atoms
	num_atoms_line = contents_dash[indices_containing_substring(contents_dash, "numAtoms")[0]].split()
	index = indices_containing_substring(contents_dash[indices_containing_substring(contents_dash, "numAtoms")[0]].split(), "numAtoms")[0]
	start = indices_containing_substring(num_atoms_line[index], '"')[0]
	end = indices_containing_substring(num_atoms_line[index], '"')[1]
	num_atoms = int(num_atoms_line[index][start+1:end])

	# Find the line where positions are given
	position_line_index = indices_containing_substring(contents_dash, "<position>")[0]

	# Find the box position
	bounds_line = contents_dash[indices_containing_substring(contents_dash, "bounds")[0]].split()

	# x lo position
	index_xlo = indices_containing_substring(bounds_line, "xlo")[0]
	start = indices_containing_substring(bounds_line[index_xlo], '"')[0]
	end = indices_containing_substring(bounds_line[index_xlo], '"')[1]
	xlo = float(bounds_line[index_xlo][start+1:end])

	# x hi position
	index_xhi = indices_containing_substring(bounds_line, "xhi")[0]
	start = indices_containing_substring(bounds_line[index_xhi], '"')[0]
	end = indices_containing_substring(bounds_line[index_xhi], '"')[1]
	xhi = float(bounds_line[index_xhi][start+1:end])

	# y lo position
	index_ylo = indices_containing_substring(bounds_line, "ylo")[0]
	start = indices_containing_substring(bounds_line[index_ylo], '"')[0]
	end = indices_containing_substring(bounds_line[index_ylo], '"')[1]
	ylo = float(bounds_line[index_ylo][start+1:end])

	# y hi position
	index_yhi = indices_containing_substring(bounds_line, "yhi")[0]
	start = indices_containing_substring(bounds_line[index_yhi], '"')[0]
	end = indices_containing_substring(bounds_line[index_yhi], '"')[1]
	yhi = float(bounds_line[index_yhi][start+1:end])

	# z lo position
	index_zlo = indices_containing_substring(bounds_line, "zlo")[0]
	start = indices_containing_substring(bounds_line[index_zlo], '"')[0]
	end = indices_containing_substring(bounds_line[index_zlo], '"')[1]
	zlo = float(bounds_line[index_zlo][start+1:end])

	# z hi position
	index_zhi = indices_containing_substring(bounds_line, "zhi")[0]
	start = indices_containing_substring(bounds_line[index_zhi], '"')[0]
	end = indices_containing_substring(bounds_line[index_zhi], '"')[1]
	zhi = float(bounds_line[index_zhi][start+1:end])

	# Calculate box dimensions
	x_boxlength = xhi - xlo
	y_boxlength = yhi - ylo
	z_boxlength = zhi - zlo

	# Atoms per molecule
	atoms_per_molecule = 20

	# List to hold positions
	positions = []

	# Compute numpy array of atomic positions
	for i in range(int(num_atoms/atoms_per_molecule)):

		molecule = []

		for j in range(atoms_per_molecule):
			molecule.append(contents_dash[position_line_index+j+1+atoms_per_molecule*i].split())

		positions.append(molecule)

	positions = np.array(positions)

	# Unwrap the molecules from periodic conditions
	for i in range(int(num_atoms/atoms_per_molecule)):

		molecule = list(positions[i, :, :])
		first_atom = molecule[0]

		for j in range(len(molecule)):

			string_contents = molecule[j]

			x_distance = float(string_contents[0]) - float(first_atom[0])
			y_distance = float(string_contents[1]) - float(first_atom[1])
			z_distance = float(string_contents[2]) - float(first_atom[2])

			if abs(x_distance) >= x_boxlength / 2.0:
				string_contents[0] = str(float(string_contents[0]) - np.sign(x_distance)*x_boxlength)

			if abs(y_distance) >= y_boxlength / 2.0:
				string_contents[1] = str(float(string_contents[1]) - np.sign(y_distance)*y_boxlength)

			if abs(z_distance) >= z_boxlength / 2.0:
				string_contents[2] = str(float(string_contents[2]) - np.sign(z_distance)*z_boxlength)

			positions[i, j, :] = string_contents


	# Update the original LAMMPS script to reflect the DASH restart coordinates
	atoms_index = indices_containing_substring(contents_lammps, "Atoms")[-1]
	atoms_list = ['' for x in range(num_atoms)]

	for i in range(num_atoms):

		string_contents = contents_lammps[atoms_index + 2 + i].split()
		string_contents[3] = str(float(string_contents[3]))
		string_contents[4] = str(float(string_contents[4]))
		string_contents[5] = str(float(string_contents[5]))
		string_contents[6] = str(float(string_contents[6]))	
		string_contents[7] = str(0)
		string_contents[8] = str(0)
		string_contents[9] = str(0)
		atoms_list[int(string_contents[0]) - 1] = '  '.join(string_contents)

	for i in range(len(atoms_list)):

		string_contents = atoms_list[i].split()
		string_contents[4] = str(float(positions[int(math.floor(i/atoms_per_molecule)), i%20, 0]))
		string_contents[5] = str(float(positions[int(math.floor(i/atoms_per_molecule)), i%20, 1]))
		string_contents[6] = str(float(positions[int(math.floor(i/atoms_per_molecule)), i%20, 2]))

		atoms_list[i] = '  '.join(string_contents)


	for p in range(num_atoms):
		contents_lammps[atoms_index + 2 + p] = '     ' + atoms_list[p] + '\n'


	# Update the box coordinates
	string_contents = contents_lammps[indices_containing_substring(contents_lammps, "xlo")[0]].split()
	string_contents[0] = str(xlo)
	string_contents[1] = str(xhi)
	contents_lammps[indices_containing_substring(contents_lammps, "xlo")[0]] = '  '.join(string_contents) + '\n'

	string_contents = contents_lammps[indices_containing_substring(contents_lammps, "ylo")[0]].split()
	string_contents[0] = str(ylo)
	string_contents[1] = str(yhi)
	contents_lammps[indices_containing_substring(contents_lammps, "ylo")[0]] = '  '.join(string_contents) + '\n'

	string_contents = contents_lammps[indices_containing_substring(contents_lammps, "zlo")[0]].split()
	string_contents[0] = str(zlo)
	string_contents[1] = str(zhi)
	contents_lammps[indices_containing_substring(contents_lammps, "zlo")[0]] = '  '.join(string_contents) + '\n'	

	# Write the final file to the output file
	out = open(files['out'], 'w')
	for k in range(len(contents_lammps)):
		out.write(contents_lammps[k])
	out.close()

		
def main(argv):

	parser = create_parser()
	args = parser.parse_args()
	files, options = convert_args(args)

	process_datafile(files)

if __name__ == "__main__":
	main(sys.argv[1:])