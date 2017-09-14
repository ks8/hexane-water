""" Script to replicate molecules for LAMMPS """
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np 
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
import math
import sys

# Open the file for writing: start with the read_data information for a single molecule of interest contained in a file named 'data.txt'
f = open('data.txt', 'r')
contents = f.readlines()
f.close()

# Set the number of molecules (must be an integer cubed) 
num_molecules = 27

# Save the original header
num_atoms_original = int(contents[2].split()[0])
num_bonds_original = int(contents[3].split()[0])
num_angles_original = int(contents[4].split()[0])
num_dihedrals_original = int(contents[5].split()[0])

# Adjust the header to reflect num_molecules
for i in range(4):
	string_contents = contents[2 + i].split()
	string_contents[0] = str(int(num_molecules*float(string_contents[0])))
	string_contents = '          ' + '  '.join(string_contents) + '\n'
	contents[2 + i] = string_contents

# Save the atomic coordinates of the original molecule
coordinates_original_matrix = np.zeros([num_atoms_original, 3])
atom_index = contents.index("Atoms\n")
for i in range(num_atoms_original):
	coordinates_contents = contents[atom_index + 2 + i].split()
	coordinates_original_matrix[i][0] = float(coordinates_contents[4])
	coordinates_original_matrix[i][1] = float(coordinates_contents[5])
	coordinates_original_matrix[i][2] = float(coordinates_contents[6])

# Calculate the center of geometry of the original molecule
center_original = np.mean(coordinates_original_matrix, axis = 0)

# Calculate the maximum dimensional range of the molecule (prevent overlapping using appropriate spacing)
ranges_original = np.max(coordinates_original_matrix, axis = 0) - np.min(coordinates_original_matrix, axis = 0)
spacing = math.ceil(np.sqrt(2.0*(max(ranges_original)**2)))

# Add new molecules to the file
atom_counter = 0
for xj in range(int(math.pow(num_molecules, 1/3))):
	for yj in range(int(math.pow(num_molecules, 1/3))):
		for zj in range(int(math.pow(num_molecules, 1/3))):
			for i in range(num_atoms_original):

				atom_counter += 1

				if(atom_counter <= num_atoms_original):
					coordinates_contents = contents[atom_index + 2 + i].split()
					coordinates_contents[4] = str(float(coordinates_contents[4]) + spacing*(float(xj) + 1.0))
					coordinates_contents[5] = str(float(coordinates_contents[5]) + spacing*(float(yj) + 1.0))
					coordinates_contents[6] = str(float(coordinates_contents[6]) + spacing*(float(zj) + 1.0))
					coordinates_contents = '     ' + '    '.join(coordinates_contents) + '\n'
					contents[atom_index + 2 + i] = coordinates_contents

				else:
					coordinates_contents = contents[atom_index + 2 + i].split()
					
					coordinates_contents[0] = str(atom_counter)
			
					coordinates_contents[4] = str(coordinates_original_matrix[i][0] + spacing*(float(xj) + 1.0))
					coordinates_contents[5] = str(coordinates_original_matrix[i][1] + spacing*(float(yj) + 1.0))
					coordinates_contents[6] = str(coordinates_original_matrix[i][2] + spacing*(float(zj) + 1.0))

					coordinates_contents = '     ' + '    '.join(coordinates_contents) + '\n'
					contents.insert(atom_index + 1 + atom_counter, coordinates_contents)

# Add new bond information to the file
bond_index = contents.index("Bonds\n")
for i in range(num_molecules - 1):

	for j in range(num_bonds_original):

		bond_contents = contents[bond_index + 2 + j].split()

		bond_contents[0] = str(int(bond_contents[0]) + num_bonds_original*(i + 1))
		bond_contents[2] = str(int(bond_contents[2]) + num_atoms_original*(i + 1))
		bond_contents[3] = str(int(bond_contents[3]) + num_atoms_original*(i + 1))

		bond_contents = '    ' + '    '.join(bond_contents) + '\n'
		contents.insert(bond_index + num_bonds_original + 2 + j + num_bonds_original*i, bond_contents)

# Add new angle information to the file
angle_index = contents.index("Angles\n")
for i in range(num_molecules - 1):
	for j in range(num_angles_original):

		angle_contents = contents[angle_index + 2 + j].split()

		angle_contents[0] = str(int(angle_contents[0]) + num_angles_original*(i + 1))

		angle_contents[2] = str(int(angle_contents[2]) + num_atoms_original*(i + 1))
		angle_contents[3] = str(int(angle_contents[3]) + num_atoms_original*(i + 1))
		angle_contents[4] = str(int(angle_contents[4]) + num_atoms_original*(i + 1))

		angle_contents = '    ' + '    '.join(angle_contents) + '\n'
		contents.insert(angle_index + num_angles_original + 2 + j + num_angles_original*i, angle_contents)

# Add new dihedral information to the file
dihedral_index = contents.index("Dihedrals\n")
for i in range(num_molecules - 1):
	for j in range(num_dihedrals_original):

		dihedral_contents = contents[dihedral_index + 2 + j].split()

		dihedral_contents[0] = str(int(dihedral_contents[0]) + num_dihedrals_original*(i + 1))

		dihedral_contents[2] = str(int(dihedral_contents[2]) + num_atoms_original*(i + 1))
		dihedral_contents[3] = str(int(dihedral_contents[3]) + num_atoms_original*(i + 1))
		dihedral_contents[4] = str(int(dihedral_contents[4]) + num_atoms_original*(i + 1))
		dihedral_contents[5] = str(int(dihedral_contents[5]) + num_atoms_original*(i + 1))

		dihedral_contents = '    ' + '    '.join(dihedral_contents) + '\n'
		contents.insert(dihedral_index + num_dihedrals_original + 2 + j + num_dihedrals_original*i, dihedral_contents)

# Rewrite file
f = open('data.txt', 'w')
contents = "".join(contents)
f.write(contents)
f.close()

# Read file
f = open('data.txt', 'r')
print(f.read())
f.close()












