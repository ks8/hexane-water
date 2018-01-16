"Script to build a full read_data LAMMPS file of the hexane water interface from two sets of molecular data files for webb hexane for the datafileoriginal input to lammps_restart.py, by Kirk Swanson"
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

# Function to parse arguments for two LAMMPS data files
def create_parser():

	# Create parser and add arguments
	parser = argparse.ArgumentParser(description='Read LAMMPS data files')
	parser.add_argument('-data1', dest='datafile1', default=None, help='Name of first LAMMPS data file')
	parser.add_argument('-N1', dest='number1', default=0, help='Number of molecules in first data file')
	parser.add_argument('-data2', dest='datafile2', default=None, help='Name of second LAMMPS data file')
	parser.add_argument('-N2', dest='number2', default=0, help='Number of molecules in second data file')
	parser.add_argument('-out', dest='outfile', default="input.txt", help='Desired output file name.  (default: input.txt')

	return parser

# Function to convert arguments into a dictionary
def convert_args(args):

	# Files dictionary
	files={}
	files['data1'] = args.datafile1
	files['data2'] = args.datafile2
	files['out'] = args.outfile

	options = {}
	options['N1'] = args.number1
	options['N2'] = args.number2

	# Print confirmation
	print("**************************")
	print("# Data file 1: ", files['data1'])
	print("# Data file 2: ", files['data2'])
	print("# Final datafile for LAMMPS input ", files['out'])
	print("**************************")
	print(" ")

	return files, options

def index_containing_substring(the_list, substring):
	for i, s in enumerate(the_list):
		if substring in s:
			return i
	return -1

def process_datafile(files, options):

	f1 = open(files['data1'], 'r')
	contents1 = f1.readlines()
	f1.close()

	f2 = open(files['data2'], 'r')
	contents2 = f2.readlines()
	f2.close()

	num_molecules_type1 = int(options['N1'])
	num_molecules_type2 = int(options['N2'])

	# Check if the second molecule contains dihedrals
	if index_containing_substring(contents2, "dihedral types") == -1:
		dihedrals = False
	else:
		dihedrals = True

	# Save the original headers for each type
	num_atoms_original_type1 = contents1[index_containing_substring(contents1, "atoms")].split()
	num_bonds_original_type1 = contents1[index_containing_substring(contents1, "bonds")].split()
	num_angles_original_type1 = contents1[index_containing_substring(contents1, "angles")].split()
	num_dihedrals_original_type1 = contents1[index_containing_substring(contents1, "dihedrals")].split()

	num_atomtype_original_type1 = contents1[index_containing_substring(contents1, "atom types")].split()
	num_bondtype_original_type1 = contents1[index_containing_substring(contents1, "bond types")].split()
	num_angletype_original_type1 = contents1[index_containing_substring(contents1, "angle types")].split()
	num_dihedraltype_original_type1 = contents1[index_containing_substring(contents1, "dihedral types")].split()

	num_atoms_original_type2 = contents2[index_containing_substring(contents2, "atoms")].split()
	num_bonds_original_type2 = contents2[index_containing_substring(contents2, "bonds")].split()
	num_angles_original_type2 = contents2[index_containing_substring(contents2, "angles")].split()
	if dihedrals == True:
		num_dihedrals_original_type2 = contents2[index_containing_substring(contents2, "dihedrals")].split()
	else:
		num_dihedrals_original_type2 = [0, 'dihedrals']

	num_atomtype_original_type2 = contents2[index_containing_substring(contents2, "atom types")].split()
	num_bondtype_original_type2 = contents2[index_containing_substring(contents2, "bond types")].split()
	num_angletype_original_type2 = contents2[index_containing_substring(contents2, "angle types")].split()
	if dihedrals == True:
		num_dihedraltype_original_type2 = contents2[index_containing_substring(contents2, "dihedral types")].split()
	else:
		num_dihedraltype_original_type2 = [0, 'dihedral types']

	# Save the atomic coordinates of the original molecule type 1
	coordinates_original_matrix_type1 = np.zeros([int(num_atoms_original_type1[0]), 3])
	atom_index_type1 = contents1.index("Atoms\n")
	for i in range(int(num_atoms_original_type1[0])):
		coordinates_contents = contents1[atom_index_type1 + 2 + i].split()
		coordinates_original_matrix_type1[i][0] = float(coordinates_contents[4])
		coordinates_original_matrix_type1[i][1] = float(coordinates_contents[5])
		coordinates_original_matrix_type1[i][2] = float(coordinates_contents[6])

	# Calculate the center of geometry of the original molecule type 1
	center_original_type1 = np.mean(coordinates_original_matrix_type1, axis = 0)

	# Calculate the maximum dimensional range of the molecule (prevent overlapping using appropriate spacing) type 1
	ranges_original_type1 = np.max(coordinates_original_matrix_type1, axis = 0) - np.min(coordinates_original_matrix_type1, axis = 0)
	spacing_type1 = math.ceil(np.sqrt(2.0*(max(ranges_original_type1)**2)))

	# Save the atomic coordinates of the original molecule type 2
	coordinates_original_matrix_type2 = np.zeros([int(num_atoms_original_type2[0]), 3])
	atom_index_type2 = contents2.index("Atoms\n")
	for i in range(int(num_atoms_original_type2[0])):
		coordinates_contents = contents2[atom_index_type2 + 2 + i].split()
		coordinates_original_matrix_type2[i][0] = float(coordinates_contents[4])
		coordinates_original_matrix_type2[i][1] = float(coordinates_contents[5])
		coordinates_original_matrix_type2[i][2] = float(coordinates_contents[6])

	# Calculate the center of geometry of the original molecule type 2
	center_original_type2 = np.mean(coordinates_original_matrix_type2, axis = 0)

	# Calculate the maximum dimensional range of the molecule (prevent overlapping using appropriate spacing) 2
	ranges_original_type2 = np.max(coordinates_original_matrix_type2, axis = 0) - np.min(coordinates_original_matrix_type2, axis = 0)
	spacing_type2 = math.ceil(np.sqrt(2.0*(max(ranges_original_type2)**2)))

	# Set the desired number of molecules for each type
	density = 1.0/pow(max(spacing_type1, spacing_type2), 3)
	volume = (1.0/density)*(num_molecules_type1 + num_molecules_type2)
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
	out.write('          ' + str(int(num_molecules_type1*float(num_atoms_original_type1[0])) + int(num_molecules_type2*float(num_atoms_original_type2[0]))) + ' atoms' + '\n')
	out.write('          ' + str(int(num_molecules_type1*float(num_bonds_original_type1[0])) + int(num_molecules_type2*float(num_bonds_original_type2[0]))) + ' bonds' + '\n')
	out.write('          ' + str(int(num_molecules_type1*float(num_angles_original_type1[0])) + int(num_molecules_type2*float(num_angles_original_type2[0]))) + ' angles' + '\n')
	out.write('          ' + str(int(num_molecules_type1*float(num_dihedrals_original_type1[0])) + int(num_molecules_type2*float(num_dihedrals_original_type2[0]))) + ' dihedrals' + '\n')
	out.write('\n')
	out.write('          ' + str(int(num_atomtype_original_type1[0]) + int(num_atomtype_original_type2[0])) + ' atom types' + '\n')
	out.write('          ' + str(int(num_bondtype_original_type1[0]) + int(num_bondtype_original_type2[0])) + ' bond types' + '\n')
	out.write('          ' + str(int(num_angletype_original_type1[0]) + int(num_angletype_original_type2[0])) + ' angle types' + '\n')
	out.write('          ' + str(int(num_dihedraltype_original_type1[0]) + int(num_dihedraltype_original_type2[0])) + ' dihedral types' + '\n')
	out.write('\n')
	out.write(str(-1.0*spacing) + ' ' + str(1.0*spacing*L) + ' ' + ' xlo xhi\n')
	out.write(str(-1.0*spacing) + ' ' + str(1.0*spacing*L) + ' ' + ' ylo yhi\n')
	out.write(str(-1.0*spacing) + ' ' + str(1.0*spacing*L) + ' ' + ' zlo zhi\n')
	out.write('\n')

	# Masses
	masses_index_type1 = contents1.index("Masses\n")
	masses_index_type2 = contents2.index("Masses\n")
	out.write('Masses\n')
	out.write('\n')
	masses_counter = 0
	for i in range(int(num_atomtype_original_type1[0])):
		masses_counter += 1
		out.write('  ' + str(masses_counter) + ' ' + str(contents1[masses_index_type1 + 2 + i].split()[1]) + '\n')
	for i in range(int(num_atomtype_original_type2[0])):
		masses_counter += 1
		out.write('  ' + str(masses_counter) + ' ' + str(contents2[masses_index_type2 + 2 + i].split()[1]) + '\n')
	out.write('\n')

	# Pair Coeffs
	# paircoeffs_index_type1 = contents1.index("Pair Coeffs\n")
	paircoeffs_index_type2 = contents2.index("Pair Coeffs\n")
	out.write('Pair Coeffs\n')
	out.write('\n')
	paircoeffs_counter = 0
	for i in range(int(num_atomtype_original_type1[0])):
		paircoeffs_counter += 1
		out.write('  ' + str(paircoeffs_counter) + '    ' + '    '.join(contents1[paircoeffs_index_type1 + 2 + i].split()[1:]) + '\n')
	for i in range(int(num_atomtype_original_type2[0])):
		paircoeffs_counter += 1
		out.write('  ' + str(paircoeffs_counter) + '    ' + '    '.join(contents2[paircoeffs_index_type2 + 2 + i].split()[1:]) + '\n')
	out.write('\n')

	# Bond Coeffs
	bondcoeffs_index_type1 = contents1.index("Bond Coeffs\n")
	bondcoeffs_index_type2 = contents2.index("Bond Coeffs\n")
	out.write('Bond Coeffs\n')
	out.write('\n')
	bondcoeffs_counter = 0
	for i in range(int(num_bondtype_original_type1[0])):
		bondcoeffs_counter += 1
		out.write('  ' + str(bondcoeffs_counter) + '    ' + '    '.join(contents1[bondcoeffs_index_type1 + 2 + i].split()[1:]) + '\n')
	for i in range(int(num_bondtype_original_type2[0])):
		bondcoeffs_counter += 1
		out.write('  ' + str(bondcoeffs_counter) + '    ' + '    '.join(contents2[bondcoeffs_index_type2 + 2 + i].split()[1:]) + '\n')
	out.write('\n')

	# Angle Coeffs
	anglecoeffs_index_type1 = contents1.index("Angle Coeffs\n")
	anglecoeffs_index_type2 = contents2.index("Angle Coeffs\n")
	out.write('Angle Coeffs\n')
	out.write('\n')
	anglecoeffs_counter = 0
	for i in range(int(num_angletype_original_type1[0])):
		anglecoeffs_counter += 1
		out.write('  ' + str(anglecoeffs_counter) + '    ' + '    '.join(contents1[anglecoeffs_index_type1 + 2 + i].split()[1:]) + '\n')
	for i in range(int(num_angletype_original_type2[0])):
		anglecoeffs_counter += 1
		out.write('  ' + str(anglecoeffs_counter) + '    ' + '    '.join(contents2[anglecoeffs_index_type2 + 2 + i].split()[1:]) + '\n')
	out.write('\n')

	# Dihedral Coeffs
	dihedralcoeffs_index_type1 = contents1.index("Dihedral Coeffs\n")
	if dihedrals == True:
		dihedralcoeffs_index_type2 = contents2.index("Dihedral Coeffs\n")
	out.write('Dihedral Coeffs\n')
	out.write('\n')
	dihedralcoeffs_counter = 0
	for i in range(int(num_dihedraltype_original_type1[0])):
		dihedralcoeffs_counter += 1
		out.write('  ' + str(dihedralcoeffs_counter) + '    ' + '    '.join(contents1[dihedralcoeffs_index_type1 + 2 + i].split()[1:]) + '\n')
	if dihedrals == True:
		for i in range(int(num_dihedraltype_original_type2[0])):
			dihedralcoeffs_counter += 1
			out.write('  ' + str(dihedralcoeffs_counter) + '    ' + '    '.join(contents2[dihedralcoeffs_index_type2 + 2 + i].split()[1:]) + '\n')
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

	# Choose molecule lattice assignments
	lattice_locations = np.random.choice(lattice, num_molecules_type1 + num_molecules_type2, replace = False)
	atom_counter = 0
	molecule_counter = 0

	for i in range(num_molecules_type1):

		lattice_point = np.array([int(x) for x in lattice_locations[molecule_counter].split()])

		molecule_counter += 1

		# Define random angles
		theta_x = np.random.uniform(0.0, math.pi)
		theta_y = np.random.uniform(0.0, math.pi)
		theta_z = np.random.uniform(0.0, math.pi)

		matrices = dict()
		matrices['x'] = rotation_matrix_x(theta_x)
		matrices['y'] = rotation_matrix_y(theta_y)
		matrices['z'] = rotation_matrix_z(theta_z)

		matrix_choices = np.random.choice(['x', 'y', 'z'], 3, replace = False)

		for j in range(int(num_atoms_original_type1[0])):

			atom_counter += 1

			coordinates_contents = contents1[atom_index_type1 + 2 + j].split()


			coordinates_contents[0] = str(atom_counter)
			coordinates_contents[1] = str(molecule_counter)

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

	for i in range(num_molecules_type2):

		lattice_point = np.array([int(x) for x in lattice_locations[molecule_counter].split()])

		molecule_counter += 1

		# Define random angles
		theta_x = np.random.uniform(0.0, math.pi)
		theta_y = np.random.uniform(0.0, math.pi)
		theta_z = np.random.uniform(0.0, math.pi)

		matrices = dict()
		matrices['x'] = rotation_matrix_x(theta_x)
		matrices['y'] = rotation_matrix_y(theta_y)
		matrices['z'] = rotation_matrix_z(theta_z)

		matrix_choices = np.random.choice(['x', 'y', 'z'], 3, replace = False)


		for j in range(int(num_atoms_original_type2[0])):

			atom_counter += 1

			coordinates_contents = contents2[atom_index_type2 + 2 + j].split()


			coordinates_contents[0] = str(atom_counter)
			coordinates_contents[1] = str(molecule_counter)
			coordinates_contents[2] = str(int(coordinates_contents[2]) + int(num_atomtype_original_type1[0]))

			new_coordinates = coordinates_original_matrix_type2[j]
			new_coordinates = new_coordinates - center_original_type2
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

	bond_index_type2 = contents2.index("Bonds\n")
	for i in range(num_molecules_type2):

		for j in range(int(num_bonds_original_type2[0])):

			bond_contents = contents2[bond_index_type2 + 2 + j].split()

			bond_contents[0] = str(int(bond_contents[0]) + int(num_bonds_original_type2[0])*(i) + num_molecules_type1*int(num_bonds_original_type1[0]))
			bond_contents[1] = str(int(bond_contents[1]) + int(num_bondtype_original_type1[0]))
			bond_contents[2] = str(int(bond_contents[2]) + int(num_atoms_original_type2[0])*(i) + num_molecules_type1*int(num_atoms_original_type1[0]))
			bond_contents[3] = str(int(bond_contents[3]) + int(num_atoms_original_type2[0])*(i) + num_molecules_type1*int(num_atoms_original_type1[0]))

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

	angle_index_type2 = contents2.index("Angles\n")
	for i in range(num_molecules_type2):

		for j in range(int(num_angles_original_type2[0])):

			angle_contents = contents2[angle_index_type2 + 2 + j].split()

			angle_contents[0] = str(int(angle_contents[0]) + int(num_angles_original_type2[0])*(i) + num_molecules_type1*int(num_angles_original_type1[0]))
			angle_contents[1] = str(int(angle_contents[1]) + int(num_angletype_original_type1[0]))
			angle_contents[2] = str(int(angle_contents[2]) + int(num_atoms_original_type2[0])*(i) + num_molecules_type1*int(num_atoms_original_type1[0]))
			angle_contents[3] = str(int(angle_contents[3]) + int(num_atoms_original_type2[0])*(i) + num_molecules_type1*int(num_atoms_original_type1[0]))
			angle_contents[4] = str(int(angle_contents[4]) + int(num_atoms_original_type2[0])*(i) + num_molecules_type1*int(num_atoms_original_type1[0]))

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

	if dihedrals == True:
		dihedral_index_type2 = contents2.index("Dihedrals\n")
		for i in range(num_molecules_type2):

			for j in range(int(num_dihedrals_original_type2[0])):

				dihedral_contents = contents2[dihedral_index_type2 + 2 + j].split()

				dihedral_contents[0] = str(int(dihedral_contents[0]) + int(num_dihedrals_original_type2[0])*(i) + num_molecules_type1*int(num_dihedrals_original_type1[0]))
				dihedral_contents[1] = str(int(dihedral_contents[1]) + int(num_dihedraltype_original_type1[0]))
				dihedral_contents[2] = str(int(dihedral_contents[2]) + int(num_atoms_original_type2[0])*(i) + num_molecules_type1*int(num_atoms_original_type1[0]))
				dihedral_contents[3] = str(int(dihedral_contents[3]) + int(num_atoms_original_type2[0])*(i) + num_molecules_type1*int(num_atoms_original_type1[0]))
				dihedral_contents[4] = str(int(dihedral_contents[4]) + int(num_atoms_original_type2[0])*(i) + num_molecules_type1*int(num_atoms_original_type1[0]))
				dihedral_contents[5] = str(int(dihedral_contents[5]) + int(num_atoms_original_type2[0])*(i) + num_molecules_type1*int(num_atoms_original_type1[0]))

				out.write('    ' + '    '.join(dihedral_contents) + '\n')


def main(argv):

	parser = create_parser()
	args = parser.parse_args()
	files, options = convert_args(args)

	process_datafile(files, options)

if __name__ == "__main__":
	main(sys.argv[1:])

