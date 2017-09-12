""" Script to replicate molecules for LAMMPS """
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np 
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Open the file for writing
f = open('data.txt', 'r')
contents = f.readlines()
atom_index = contents.index("Atoms\n")
f.close()

# Header
for i in range(4):
	string_contents = contents[2 + i].split()
	string_contents[0] = str(int(2.0*float(string_contents[0])))
	string_contents = '          ' + '  '.join(string_contents) + '\n'
	contents[2 + i] = string_contents

# Atom positions
for i in range(20):
	string_contents = contents[atom_index + 2 + i].split()
	string_contents[0] = str(int(string_contents[0]) + 20)
	string_contents[4] = str(float(string_contents[4]) + 5.00000)
	string_contents[5] = str(float(string_contents[5]) + 5.00000)
	string_contents[6] = str(float(string_contents[6]) + 5.00000)
	string_contents = '     ' + '    '.join(string_contents) + '\n'
	contents.insert(atom_index + 22 + i, string_contents)

bond_index = contents.index("Bonds\n")
for i in range(19):
	string_contents = contents[bond_index + 2 + i].split()
	string_contents[0] = str(int(string_contents[0]) + 19)
	string_contents[2] = str(int(string_contents[2]) + 20)
	string_contents[3] = str(int(string_contents[3]) + 20)
	string_contents = '    ' + '    '.join(string_contents) + '\n'
	contents.insert(bond_index + 21 + i, string_contents)

angle_index = contents.index("Angles\n")
for i in range(36):
	string_contents = contents[angle_index + 2 + i].split()
	string_contents[0] = str(int(string_contents[0]) + 36)
	string_contents[2] = str(int(string_contents[2]) + 20)
	string_contents[3] = str(int(string_contents[3]) + 20)
	string_contents[4] = str(int(string_contents[4]) + 20)
	string_contents = '    ' + '    '.join(string_contents) + '\n'
	contents.insert(angle_index + 38 + i, string_contents)

dihedral_index = contents.index("Dihedrals\n")
for i in range(45):
	string_contents = contents[dihedral_index + 2 + i].split()
	string_contents[0] = str(int(string_contents[0]) + 45)
	string_contents[2] = str(int(string_contents[2]) + 20)
	string_contents[3] = str(int(string_contents[3]) + 20)
	string_contents[4] = str(int(string_contents[4]) + 20)
	string_contents[5] = str(int(string_contents[5]) + 20)
	string_contents = '    ' + '    '.join(string_contents) + '\n'
	contents.insert(dihedral_index + 47 + i, string_contents)

# Rewrite file
f = open('data.txt', 'w')
contents = "".join(contents)
f.write(contents)
f.close()

# Read file
f = open('data.txt', 'r')
print(f.read())
f.close()












