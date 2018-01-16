"Script to build a full read_data LAMMPS file from previous simulation outputs, by Kirk Swanson"
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

# Function to parse arguments for one LAMMPS data files
def create_parser():

	# Create parser and add arguments
	parser = argparse.ArgumentParser(description='Read data files')
	parser.add_argument('-data1', dest='datafile1', default=None, help='Name of first LAMMPS output data file')
	parser.add_argument('-data2', dest='datafile2', default=None, help='Name of second LAMMPS output data file')
	parser.add_argument('-dataoriginal', dest='datafileoriginal', default=None, help='Name of original LAMMPS input data file')
	parser.add_argument('-out', dest='outfile', default="input_restart.txt", help='Desired output file name.  (default: input_restart.txt')
	
	return parser

# Function to convert arguments into a dictionary
def convert_args(args):

	# Files dictionary
	files={}
	files['data1'] = args.datafile1
	files['data2'] = args.datafile2
	files['dataoriginal'] = args.datafileoriginal
	files['out'] = args.outfile
	
	options = {}

	# Print confirmation
	print("**************************")
	print("# Data file 1: ", files['data1'])
	print("# Data file 2: ", files['data2'])
	print("# Data file original: ", files['dataoriginal'])
	print("# Final datafile for LAMMPS input ", files['out'])
	print("**************************")
	print(" ")

	return files, options

def indices_containing_substring(the_list, substring):
	indices = []
	for i, s in enumerate(the_list):
		if substring in s:
			indices.append(i)	
	return indices

def process_datafile(files):

	f1 = open(files['data1'], 'r')
	contents1 = f1.readlines()
	f1.close()

	f2 = open(files['data2'], 'r')
	contents2 = f2.readlines()
	f2.close()

	f_original = open(files['dataoriginal'], 'r')
	contents_original = f_original.readlines()
	f_original.close()
	
	num_atoms_type1 = int(contents1[indices_containing_substring(contents1, "atoms")[0]].split()[0]) 	
	atoms_index_type1 = indices_containing_substring(contents1, "Atoms")[-1]
	atoms_type1 = ['' for x in range(num_atoms_type1)]

	for i in range(num_atoms_type1):

		string_contents = contents1[atoms_index_type1 + 2 + i].split()
		string_contents[3] = str(float(string_contents[3]))
		string_contents[4] = str(float(string_contents[4]))
		string_contents[5] = str(float(string_contents[5]))
		string_contents[6] = str(float(string_contents[6]))	
		string_contents[7] = str(0)
		string_contents[8] = str(0)
		string_contents[9] = str(0)
		atoms_type1[int(string_contents[0]) - 1] = '  '.join(string_contents)

	num_molecules_type1 = int(atoms_type1[num_atoms_type1 - 1].split()[1])

	x_boxlength = float(contents1[indices_containing_substring(contents1, "xlo")[0]].split()[1]) - float(contents1[indices_containing_substring(contents1, "xlo")[0]].split()[0])
	y_boxlength = float(contents1[indices_containing_substring(contents1, "ylo")[0]].split()[1]) - float(contents1[indices_containing_substring(contents1, "ylo")[0]].split()[0])
	z_boxlength = float(contents1[indices_containing_substring(contents1, "zlo")[0]].split()[1]) - float(contents1[indices_containing_substring(contents1, "zlo")[0]].split()[0])

	for i in range(num_molecules_type1):
	 	
	  	molecule = [row for row in atoms_type1 if int(row.split()[1]) == i+1]
	  	first_atom = molecule[0]
	  	for j in range(len(molecule)):
	  		x_distance = float(molecule[j].split()[4]) - float(first_atom.split()[4])
	  		y_distance = float(molecule[j].split()[5]) - float(first_atom.split()[5])
	  		z_distance = float(molecule[j].split()[6]) - float(first_atom.split()[6])

	  		string_contents = molecule[j].split()
	  		
	  		if abs(x_distance) >= x_boxlength / 2.0:
				string_contents[4] = str(float(string_contents[4]) - np.sign(x_distance)*x_boxlength)

			if abs(y_distance) >= y_boxlength / 2.0:
				string_contents[5] = str(float(string_contents[5]) - np.sign(y_distance)*y_boxlength)

			if abs(z_distance) >= z_boxlength / 2.0:
				string_contents[6] = str(float(string_contents[6]) - np.sign(z_distance)*z_boxlength)

	  		atoms_type1[int(string_contents[0]) - 1] = '  '.join(string_contents)

		
	num_atoms_type2 = int(contents2[indices_containing_substring(contents2, "atoms")[0]].split()[0]) 	
	atoms_index_type2 = indices_containing_substring(contents2, "Atoms")[-1]
	atoms_type2 = ['' for x in range(num_atoms_type2)]

	for i in range(num_atoms_type2):
		string_contents = contents2[atoms_index_type2 + 2 + i].split()
		string_contents[0] = str(int(string_contents[0]) + num_atoms_type1)
		string_contents[1] = str(int(string_contents[1]) + 500)
		string_contents[2] = str(int(string_contents[2]) + 5)
		string_contents[3] = str(float(string_contents[3]))
		string_contents[4] = str(float(string_contents[4]) + 60.808 + 5)
		string_contents[5] = str(float(string_contents[5]) + 12.478)
		string_contents[6] = str(float(string_contents[6]) + 12.478)	
		string_contents[7] = str(0)
		string_contents[8] = str(0)
		string_contents[9] = str(0)
		atoms_type2[int(string_contents[0]) - 1 - num_atoms_type1] = '  '.join(string_contents)

	num_molecules_type2 = int(atoms_type2[num_atoms_type2 - 1].split()[1]) - num_molecules_type1

	x_boxlength = float(contents2[indices_containing_substring(contents2, "xlo")[0]].split()[1]) - float(contents2[indices_containing_substring(contents2, "xlo")[0]].split()[0])
	y_boxlength = float(contents2[indices_containing_substring(contents2, "ylo")[0]].split()[1]) - float(contents2[indices_containing_substring(contents2, "ylo")[0]].split()[0])
	z_boxlength = float(contents2[indices_containing_substring(contents2, "zlo")[0]].split()[1]) - float(contents2[indices_containing_substring(contents2, "zlo")[0]].split()[0])


	for i in range(num_molecules_type2):
	 	
	  	molecule = [row for row in atoms_type2 if int(row.split()[1]) == i+1+500]
	  	first_atom = molecule[0]
	  	for j in range(len(molecule)):
	  		x_distance = float(molecule[j].split()[4]) - float(first_atom.split()[4])
	  		y_distance = float(molecule[j].split()[5]) - float(first_atom.split()[5])
	  		z_distance = float(molecule[j].split()[6]) - float(first_atom.split()[6])

	  		string_contents = molecule[j].split()
	  		
	  		if abs(x_distance) >= x_boxlength / 2.0:
				string_contents[4] = str(float(string_contents[4]) - np.sign(x_distance)*x_boxlength)

			if abs(y_distance) >= y_boxlength / 2.0:
				string_contents[5] = str(float(string_contents[5]) - np.sign(y_distance)*y_boxlength)

			if abs(z_distance) >= z_boxlength / 2.0:
				string_contents[6] = str(float(string_contents[6]) - np.sign(z_distance)*z_boxlength)

	  		atoms_type2[int(string_contents[0]) - 1 - num_atoms_type1] = '  '.join(string_contents)


	atoms_index_original = indices_containing_substring(contents_original, "Atoms")[-1]

	for p in range(num_atoms_type1):
		contents_original[atoms_index_original + 2 + p] = '     ' + atoms_type1[p] + '\n'

	for p in range(num_atoms_type2):
		contents_original[atoms_index_original + 2 + num_atoms_type1 + p] = '     ' + atoms_type2[p] + '\n'

	xlo = float(contents_original[atoms_index_original + 2].split()[4])
	xhi = float(contents_original[atoms_index_original + 2].split()[4]) 
	ylo = float(contents_original[atoms_index_original + 2].split()[5]) 
	yhi = float(contents_original[atoms_index_original + 2].split()[5]) 
	zlo = float(contents_original[atoms_index_original + 2].split()[6]) 
	zhi = float(contents_original[atoms_index_original + 2].split()[6]) 

	for i in range(num_atoms_type1 + num_atoms_type2):
		x_pos = float(contents_original[atoms_index_original + 2 + i].split()[4])
		y_pos = float(contents_original[atoms_index_original + 2 + i].split()[5])
		z_pos = float(contents_original[atoms_index_original + 2 + i].split()[6])

		if x_pos > xhi:
			xhi = x_pos
		if x_pos < xlo:
			xlo = x_pos
		
		if y_pos > yhi:
			yhi = y_pos
		if y_pos < ylo:
			ylo = y_pos

		if z_pos > zhi:
			zhi = z_pos
		if z_pos < zlo:
			zlo = z_pos

	string_contents = contents_original[indices_containing_substring(contents_original, "xlo")[0]].split()
	string_contents[0] = str(xlo)
	string_contents[1] = str(xhi)
	contents_original[indices_containing_substring(contents_original, "xlo")[0]] = ' '.join(string_contents) + '\n'

	string_contents = contents_original[indices_containing_substring(contents_original, "ylo")[0]].split()
	string_contents[0] = str(ylo)
	string_contents[1] = str(yhi)
	contents_original[indices_containing_substring(contents_original, "ylo")[0]] = ' '.join(string_contents) + '\n'

	string_contents = contents_original[indices_containing_substring(contents_original, "zlo")[0]].split()
	string_contents[0] = str(zlo)
	string_contents[1] = str(zhi)
	contents_original[indices_containing_substring(contents_original, "zlo")[0]] = ' '.join(string_contents) + '\n'	





	out = open(files['out'], 'w')
	for k in range(len(contents_original)):
		out.write(contents_original[k])
	out.close()

		
def main(argv):

	parser = create_parser()
	args = parser.parse_args()
	files, options = convert_args(args)

	process_datafile(files)

if __name__ == "__main__":
	main(sys.argv[1:])