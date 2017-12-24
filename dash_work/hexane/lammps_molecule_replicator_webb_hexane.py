"Script to build a full read_data LAMMPS file from a LAMMPS data file of a single molecule of Webb hexane, by Kirk Swanson"
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

# Function to parse arguments 
def create_parser():

	# Create parser and add arguments
	parser = argparse.ArgumentParser(description='Read LAMMPS data files')
	parser.add_argument('-datafile', dest='datafile', default=None, help='Name of first LAMMPS data file')
	parser.add_argument('-num_molecules', dest='num_molecules', default=None, help='Desired number of molecules')
	parser.add_argument('-out', dest='outfile', default="input.txt", help='Desired output file name.  (default: input.txt')

	return parser

# Function to convert arguments into a dictionary
def convert_args(args):

	# Files dictionary
	files={}
	files['datafile'] = args.datafile
	files['out'] = args.outfile

	options = {}
	options['num_molecules'] = args.num_molecules

	# Print confirmation
	print("**************************")
	print("# Data file 1: ", files['datafile'])
	print("# Final datafile for LAMMPS input ", files['out'])
	print("**************************")
	print(" ")

	return files, options

def process_datafile(files, options):

	f1 = open(files['datafile'], 'r')
	contents1 = f1.readlines()
	f1.close()

	# Save the original header
	num_atoms_original_type1 = contents1[2].split()
	num_bonds_original_type1 = contents1[3].split()
	num_angles_original_type1 = contents1[4].split()
	num_dihedrals_original_type1 = contents1[5].split()

	num_atomtype_original_type1 = contents1[7].split()
	num_bondtype_original_type1 = contents1[8].split()
	num_angletype_original_type1 = contents1[9].split()
	num_dihedraltype_original_type1 = contents1[10].split()

	# Save the atomic coordinates of the original molecule 
	coordinates_original_matrix_type1 = np.zeros([int(num_atoms_original_type1[0]), 3])
	atom_index_type1 = contents1.index("Atoms\n")
	for i in range(int(num_atoms_original_type1[0])):
		coordinates_contents = contents1[atom_index_type1 + 2 + i].split()
		coordinates_original_matrix_type1[i][0] = float(coordinates_contents[4])
		coordinates_original_matrix_type1[i][1] = float(coordinates_contents[5])
		coordinates_original_matrix_type1[i][2] = float(coordinates_contents[6])

	# Calculate the center of geometry of the original molecule 
	center_original_type1 = np.mean(coordinates_original_matrix_type1, axis = 0)

	# Calculate the maximum dimensional range of the molecule (prevent overlapping using appropriate spacing) 
	ranges_original_type1 = np.max(coordinates_original_matrix_type1, axis = 0) - np.min(coordinates_original_matrix_type1, axis = 0)
	spacing_type1 = math.ceil(np.sqrt(2.0*(max(ranges_original_type1)**2)))

	# Set the desired number of molecules 
	density = 1.0/pow(spacing_type1, 3)
	num_molecules_type1 = int(options['num_molecules'])
	volume = (1.0/density)*(num_molecules_type1)
	spacing = int(math.ceil(pow(1.0/density, 1/3)))
	L = int(math.ceil(int(math.ceil(pow(volume, 1/3)))/spacing))

	# Construct a set of lattice points corresponding to the desired density and number of molecules
	lattice = []
	for i in range(L):
		for j in range(L):
			for k in range(L):
				lattice.append(str(spacing*i) + ' ' + str(spacing*j) + ' ' + str(spacing*k))	


	# Construct the header
	out = open(files['out'], 'w')
	out.write('Multiple Molecule Type Simulation\n')
	out.write('\n')
	out.write('          ' + str(int(num_molecules_type1*float(num_atoms_original_type1[0]))) + ' atoms' + '\n')
	out.write('          ' + str(int(num_molecules_type1*float(num_bonds_original_type1[0]))) + ' bonds' + '\n')
	out.write('          ' + str(int(num_molecules_type1*float(num_angles_original_type1[0]))) + ' angles' + '\n')
	out.write('          ' + str(int(num_molecules_type1*float(num_dihedrals_original_type1[0]))) + ' dihedrals' + '\n')
	out.write('\n')
	out.write('          ' + str(int(num_atomtype_original_type1[0])) + ' atom types' + '\n')
	out.write('          ' + str(int(num_bondtype_original_type1[0])) + ' bond types' + '\n')
	out.write('          ' + str(int(num_angletype_original_type1[0])) + ' angle types' + '\n')
	out.write('          ' + str(int(num_dihedraltype_original_type1[0])) + ' dihedral types' + '\n')
	out.write('\n')
	out.write(str(-1.0*spacing) + ' ' + str(1.0*spacing*L) + ' ' + ' xlo xhi\n')
	out.write(str(-1.0*spacing) + ' ' + str(1.0*spacing*L) + ' ' + ' ylo yhi\n')
	out.write(str(-1.0*spacing) + ' ' + str(1.0*spacing*L) + ' ' + ' zlo zhi\n')
	out.write('\n')

	# Write the Masses
	masses_index_type1 = contents1.index("Masses\n")
	out.write('Masses\n')
	out.write('\n')
	masses_counter = 0
	for i in range(int(num_atomtype_original_type1[0])):
		masses_counter += 1
		out.write('  ' + str(masses_counter) + ' ' + str(contents1[masses_index_type1 + 2 + i].split()[1]) + '\n')
	out.write('\n')


	# Define rotation matrices
	def rotation_matrix_x(theta):
		rotation_matrix = np.zeros([3,3])
		rotation_matrix[0][0] = 1.0
		rotation_matrix[0][1] = 0.0
		rotation_matrix[0][2] = 0.0
		rotation_matrix[1][0] = 0.0
		rotation_matrix[1][1] = math.cos(theta)
		rotation_matrix[1][2] = -1.0*math.sin(theta)
		rotation_matrix[2][0] = 0.0
		rotation_matrix[2][1] = math.sin(theta)
		rotation_matrix[2][2] = math.cos(theta)

		return rotation_matrix

	def rotation_matrix_y(theta):
		rotation_matrix = np.zeros([3,3])
		rotation_matrix[0][0] = math.cos(theta)
		rotation_matrix[0][1] = 0.0
		rotation_matrix[0][2] = math.sin(theta)
		rotation_matrix[1][0] = 0.0
		rotation_matrix[1][1] = 1.0
		rotation_matrix[1][2] = 0.0
		rotation_matrix[2][0] = -1.0*math.sin(theta)
		rotation_matrix[2][1] = 0.0
		rotation_matrix[2][2] = math.cos(theta)

		return rotation_matrix

	def rotation_matrix_z(theta):
		rotation_matrix = np.zeros([3,3])
		rotation_matrix[0][0] = math.cos(theta)
		rotation_matrix[0][1] = -1.0*math.sin(theta)
		rotation_matrix[0][2] = 0.0
		rotation_matrix[1][0] = math.sin(theta)
		rotation_matrix[1][1] = math.cos(theta)
		rotation_matrix[1][2] = 0.0
		rotation_matrix[2][0] = 0.0
		rotation_matrix[2][1] = 0.0
		rotation_matrix[2][2] = 1.0

		return rotation_matrix

	# Atoms
	out.write('Atoms\n')
	out.write('\n')

	# Choose molecule lattice assignments and set up an atom and molecule counter
	lattice_locations = np.random.choice(lattice, num_molecules_type1, replace = False)
	atom_counter = 0
	molecule_counter = 0

	# Set molecules on the lattice
	for i in range(num_molecules_type1):

		# Choose a lattice point
		lattice_point = np.array([int(x) for x in lattice_locations[molecule_counter].split()])

		# Update the molecule counter
		molecule_counter += 1

		# Define random angles
		theta_x = np.random.uniform(0.0, math.pi)
		theta_y = np.random.uniform(0.0, math.pi)
		theta_z = np.random.uniform(0.0, math.pi)

		# Construct a rotation matrix
		matrices = dict()
		matrices['x'] = rotation_matrix_x(theta_x)
		matrices['y'] = rotation_matrix_y(theta_y)
		matrices['z'] = rotation_matrix_z(theta_z)

		# Specify the rotation matrix
		matrix_choices = np.random.choice(['x', 'y', 'z'], 3, replace = False)

		# Set the atoms of this molecule at a point on the lattice with appropriate random rotation
		for j in range(int(num_atoms_original_type1[0])):

			# Update the atom counter
			atom_counter += 1

			# Select the original atom coordinates from the original data file
			coordinates_contents = contents1[atom_index_type1 + 2 + j].split()

			# Update the atom and molecule number
			coordinates_contents[0] = str(atom_counter)
			coordinates_contents[1] = str(molecule_counter)

			# Update coordinates for the lattice position with rotation
			new_coordinates = coordinates_original_matrix_type1[j]
			new_coordinates = new_coordinates - center_original_type1
			new_coordinates = np.dot(matrices[matrix_choices[0]], new_coordinates)
			new_coordinates = np.dot(matrices[matrix_choices[1]], new_coordinates)
			new_coordinates = np.dot(matrices[matrix_choices[2]], new_coordinates)
			new_coordinates = new_coordinates + lattice_point

			coordinates_contents[4] = str(float(new_coordinates[0]))
			coordinates_contents[5] = str(float(new_coordinates[1]))
			coordinates_contents[6] = str(float(new_coordinates[2]))

			out.write('     ' + '    '.join(coordinates_contents) + '\n')

	# Bonds
	out.write('\n')
	out.write('Bonds\n')
	out.write('\n')

	# Add new bond information to the file
	bond_index_type1 = contents1.index("Bonds\n")
	for i in range(num_molecules_type1):

		for j in range(int(num_bonds_original_type1[0])):

			bond_contents = contents1[bond_index_type1 + 2 + j].split()

			bond_contents[0] = str(int(bond_contents[0]) + int(num_bonds_original_type1[0])*(i))
			bond_contents[2] = str(int(bond_contents[2]) + int(num_atoms_original_type1[0])*(i))
			bond_contents[3] = str(int(bond_contents[3]) + int(num_atoms_original_type1[0])*(i))

			out.write('    ' + '    '.join(bond_contents) + '\n')

	# Angles
	out.write('\n')
	out.write('Angles\n')
	out.write('\n')

	# Add new angle information to the file
	angle_index_type1 = contents1.index("Angles\n")
	for i in range(num_molecules_type1):

		for j in range(int(num_angles_original_type1[0])):

			angle_contents = contents1[angle_index_type1 + 2 + j].split()

			angle_contents[0] = str(int(angle_contents[0]) + int(num_angles_original_type1[0])*(i))
			angle_contents[2] = str(int(angle_contents[2]) + int(num_atoms_original_type1[0])*(i))
			angle_contents[3] = str(int(angle_contents[3]) + int(num_atoms_original_type1[0])*(i))
			angle_contents[4] = str(int(angle_contents[4]) + int(num_atoms_original_type1[0])*(i))

			out.write('    ' + '    '.join(angle_contents) + '\n')

	# Dihedrals
	out.write('\n')
	out.write('Dihedrals\n')
	out.write('\n')

	# Add new dihedral information to the file
	dihedral_index_type1 = contents1.index("Dihedrals\n")
	for i in range(num_molecules_type1):

		for j in range(int(num_dihedrals_original_type1[0])):

			dihedral_contents = contents1[dihedral_index_type1 + 2 + j].split()

			dihedral_contents[0] = str(int(dihedral_contents[0]) + int(num_dihedrals_original_type1[0])*(i))
			dihedral_contents[2] = str(int(dihedral_contents[2]) + int(num_atoms_original_type1[0])*(i))
			dihedral_contents[3] = str(int(dihedral_contents[3]) + int(num_atoms_original_type1[0])*(i))
			dihedral_contents[4] = str(int(dihedral_contents[4]) + int(num_atoms_original_type1[0])*(i))
			dihedral_contents[5] = str(int(dihedral_contents[5]) + int(num_atoms_original_type1[0])*(i))

			out.write('    ' + '    '.join(dihedral_contents) + '\n')


def main(argv):

	parser = create_parser()
	args = parser.parse_args()
	files, options = convert_args(args)

	process_datafile(files, options)

if __name__ == "__main__":
	main(sys.argv[1:])

