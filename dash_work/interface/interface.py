"Script to perform a DASH hexane-water interface simulation from a DASH restart file of TIP4P/F water and a LAMMPS input file of Webb hexane, by Kirk Swanson"
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

sys.path = sys.path + ['/home/swansonk1/OLD_DASH/md_engine/build/python/build/lib.linux-x86_64-2.7']
sys.path.append('/home/swansonk1/OLD_DASH/md_engine/util_py')

from DASH import *
from LAMMPS_Reader import LAMMPS_Reader
import re
import water
from water import *

# Function to parse arguments for one DASH restart file
def create_parser():

    # Create parser and add arguments
    parser = argparse.ArgumentParser(description='Read data files')
    parser.add_argument('-water_restart_file', dest='water_restart_file', default=None, help='Name of water DASH restart xml file')
    parser.add_argument('-hexane_lammps_file', dest='hexane_lammps_file', default=None, help='Name of hexane LAMMPS input file')

    return parser

# Function to convert arguments into a dictionary
def convert_args(args):

    # Files dictionary
    files={}
    files['water_restart_file'] = args.water_restart_file
    files['hexane_lammps_file'] = args.hexane_lammps_file

    options = {}

    # Print confirmation
    print("**************************")
    print("# Data file: ", files['water_restart_file'])
    print("# Data file: ", files['hexane_lammps_file'])
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


	# Read restart files
	f1 = open(files['water_restart_file'], 'r')
	contents1 = f1.readlines()
	f1.close()

	f2 = open(files['hexane_lammps_file'], 'r')
	contents_hexane = f2.readlines()
	f2.close()

	# Extract the box coordinates from the hexane file
	hexane_x_bounds = contents_hexane[indices_containing_substring(contents_hexane, "xlo")[0]].split()
	hexane_xlo = float(hexane_x_bounds[0])
	hexane_xhi = float(hexane_x_bounds[1])

	hexane_y_bounds = contents_hexane[indices_containing_substring(contents_hexane, "ylo")[0]].split()
	hexane_ylo = float(hexane_y_bounds[0])
	hexane_yhi = float(hexane_y_bounds[1])

	hexane_z_bounds = contents_hexane[indices_containing_substring(contents_hexane, "zlo")[0]].split()
	hexane_zlo = float(hexane_z_bounds[0])
	hexane_zhi = float(hexane_z_bounds[1])

	# Save the  header of the hexane file
	num_atoms_hexane = contents_hexane[2].split()
	num_bonds_hexane = contents_hexane[3].split()
	num_angles_hexane = contents_hexane[4].split()
	num_dihedrals_hexane = contents_hexane[5].split()

	num_atomtype_hexane = contents_hexane[7].split()
	num_bondtype_hexane = contents_hexane[8].split()
	num_angletype_hexane = contents_hexane[9].split()
	num_dihedraltype_hexane = contents_hexane[10].split()

	# Save the atomic coordinates of the original molecule
	coordinates_hexane = np.zeros([int(num_atoms_hexane[0]), 3])
	atom_index_hexane = contents_hexane.index("Atoms\n")
	for i in range(int(num_atoms_hexane[0])):
		coordinates_contents = contents_hexane[atom_index_hexane + 2 + i].split()
		coordinates_hexane[i][0] = float(coordinates_contents[4])
		coordinates_hexane[i][1] = float(coordinates_contents[5])
		coordinates_hexane[i][2] = float(coordinates_contents[6])

	# Calculate the center of geometry of the original molecule
	center_original_type1 = np.mean(coordinates_hexane, axis = 0)

	# Atoms per molecule for hexane
	atoms_per_molecule_hexane = 20

	# Calculate the maximum dimensional range of the molecule (prevent overlapping using appropriate spacing)
	ranges_hexane = max(np.max(coordinates_hexane[0:atoms_per_molecule_hexane], axis = 0) - np.min(coordinates_hexane[0:atoms_per_molecule_hexane], axis = 0))


	# Calculate number of atoms
	num_atoms_line = contents1[indices_containing_substring(contents1, "numAtoms")[0]].split()
	index = indices_containing_substring(contents1[indices_containing_substring(contents1, "numAtoms")[0]].split(), "numAtoms")[0]
	start = indices_containing_substring(num_atoms_line[index], '"')[0]
	end = indices_containing_substring(num_atoms_line[index], '"')[1]
	num_atoms = int(num_atoms_line[index][start+1:end])

	# Find the line where positions are given
	position_line_index = indices_containing_substring(contents1, "<position>")[0]

	# Atoms per molecule
	atoms_per_molecule = 4

	# Find the box position
	bounds_line = contents1[indices_containing_substring(contents1, "bounds")[0]].split()

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

	x_boxlength = xhi - xlo
	y_boxlength = yhi - ylo
	z_boxlength = zhi - zlo


	# List to hold positions
	positions = []

	# List to hold atom positions
	for i in range(int(num_atoms/atoms_per_molecule)):

		molecule = []
		molecule.append(contents1[position_line_index+1+atoms_per_molecule*i].split())
		molecule.append(contents1[position_line_index+2+atoms_per_molecule*i].split())
		molecule.append(contents1[position_line_index+3+atoms_per_molecule*i].split())
		molecule.append(contents1[position_line_index+4+atoms_per_molecule*i].split())

		positions.append(molecule)

	positions = np.array(positions)

	# Move the box
	for i in range(int(num_atoms/atoms_per_molecule)):
		for j in range(atoms_per_molecule):
			positions[i,j,0] = str(float(positions[i,j,0]) + hexane_xlo - xlo)
			positions[i,j,1] = str(float(positions[i,j,1]) + hexane_ylo - ylo)
			positions[i,j,2] = str(float(positions[i,j,2]) - (zhi - hexane_zlo + ranges_hexane))

	# Unwrap molecules
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



	state = State()
	state.deviceManager.setDevice(0)
	state.rCut = 12.0
	state.padding = 2.0
	state.periodicInterval = 1
	state.shoutEvery = 100
	state.units.setReal()
	state.dt = 0.5
	nSteps   = 400000
	the_temp = 298.15
	printFreq= 1000

	oxygenHandle = 'OW'
	hydrogenHandle = 'HY'
	mSiteHandle = 'M'

	# nBeads = 1
	# the_temp = 298.15
	# the_temp *= nBeads

	state.atomParams.addSpecies(handle=oxygenHandle, mass=15.9994, atomicNum=8)
	state.atomParams.addSpecies(handle=hydrogenHandle, mass=1.0074, atomicNum=1)
	state.atomParams.addSpecies(handle=mSiteHandle,mass=0.0,atomicNum=0)

	ljcut = FixLJCut(state, 'ljcut')
	bondHARM = FixBondHarmonic(state, 'bondHARM')
	angleHARM  = FixAngleHarmonic(state, 'angleHARM')
	dihedralOPLS = FixDihedralOPLS(state, 'dihedralOPLS')

	epsilon = 0.1852 # given in kcal/mol
	sigma = 3.1589 # given in Angstroms

	ljcut.setParameter('sig',oxygenHandle, oxygenHandle, sigma)
	ljcut.setParameter('eps',oxygenHandle, oxygenHandle, epsilon)
	ljcut.setParameter('sig',hydrogenHandle, hydrogenHandle, 0.0)
	ljcut.setParameter('eps',hydrogenHandle, hydrogenHandle, 0.0)
	ljcut.setParameter('sig',mSiteHandle, mSiteHandle, 0.0)
	ljcut.setParameter('eps',mSiteHandle, mSiteHandle, 0.0)

	flexibleTIP4P = FixTIP4PFlexible(state,'TIP4PFlexible')

	bondQuart = FixBondQuartic(state,'bondQuart')
	bondQuart.setBondTypeCoefs(type=0,k2=607.19,k3=-1388.65,k4=1852.58,r0=0.9419)

	angleHARM.setAngleTypeCoefs(type=0, k=87.85, theta0=( (107.4/180.0) * pi))

	positions = positions.astype(float)

	for i in range(3650):
		center = positions[i, :, :]

		molecule = create_TIP4P_Flexible(state,oxygenHandle,hydrogenHandle,mSiteHandle,center,"restart")

		ids = []
		for atomId in molecule.ids:
			ids.append(atomId)

		# sanity check: set the mass of the M-site to zero, again
		state.atoms[ids[-1]].mass = 0.00
		# add the atom ids to our instance of FixTIP4PFlexible
		flexibleTIP4P.addMolecule(ids[0], ids[1], ids[2], ids[3])

		bondQuart.createBond(state.atoms[ids[0]], state.atoms[ids[1]],type=0)
		bondQuart.createBond(state.atoms[ids[0]], state.atoms[ids[2]],type=0)

		angleHARM.createAngle(state.atoms[ids[1]], state.atoms[ids[0]], state.atoms[ids[2]],type=0)


	state.activateFix(ljcut)
	state.activateFix(bondHARM)
	state.activateFix(angleHARM)
	state.activateFix(dihedralOPLS)
	state.activateFix(flexibleTIP4P)
	state.activateFix(bondQuart)


	charge = FixChargeEwald(state, 'charge', 'all')
	charge.setError(0.01,state.rCut,3)
	state.activateFix(charge)

	fixNVT = FixNVTAndersen(state, handle='nvt', groupHandle='all', temp=the_temp, nu=0.01)
	fixNVT.setParameters(10)
	state.activateFix(fixNVT)

	fixPressure = FixPressureBerendsen(state, 'npt', 1.0, 1000, 5)
	state.activateFix(fixPressure)



	reader = LAMMPS_Reader(state=state, setBounds=True, nonbondFix = ljcut, bondFix = bondHARM, angleFix = angleHARM, dihedralFix = dihedralOPLS, atomTypePrefix = 'lmps_')
	reader.read(dataFn = 'hexane_restart.txt', inputFns = ['hexane.in.settings'])


	# Establish geometric mixing rules
	epsilon_lmps_0 = 0.018252
	sigma_lmps_0 = 2.467643

	epsilon_lmps_1 = 0.010204
	sigma_lmps_1 = 2.447434

	epsilon_lmps_2 = 0.129282
	sigma_lmps_2 = 3.386973

	epsilon_lmps_3 = 0.133852
	sigma_lmps_3 = 3.392671

	epsilon_lmps_4 = 0.130579
	sigma_lmps_4 = 3.430693


	ljcut.setParameter('eps',oxygenHandle, 'lmps_0', math.sqrt(epsilon*epsilon_lmps_0))
	ljcut.setParameter('sig',oxygenHandle, 'lmps_0', (sigma + sigma_lmps_0)/2)

	ljcut.setParameter('eps',oxygenHandle, 'lmps_1', math.sqrt(epsilon*epsilon_lmps_1))
	ljcut.setParameter('sig',oxygenHandle, 'lmps_1', (sigma + sigma_lmps_1)/2)

	ljcut.setParameter('eps',oxygenHandle, 'lmps_2', math.sqrt(epsilon*epsilon_lmps_2))
	ljcut.setParameter('sig',oxygenHandle, 'lmps_2', (sigma + sigma_lmps_2)/2)

	ljcut.setParameter('eps',oxygenHandle, 'lmps_3', math.sqrt(epsilon*epsilon_lmps_3))
	ljcut.setParameter('sig',oxygenHandle, 'lmps_3', (sigma + sigma_lmps_3)/2)

	ljcut.setParameter('eps',oxygenHandle, 'lmps_4', math.sqrt(epsilon*epsilon_lmps_4))
	ljcut.setParameter('sig',oxygenHandle, 'lmps_4', (sigma + sigma_lmps_4)/2)


	# ljcut.setParameter('eps',oxygenHandle, 'lmps_0', 2*math.sqrt(epsilon*epsilon_lmps_0)*((pow(sigma, 3) + pow(sigma_lmps_0, 3))/(pow(sigma, 6) + pow(sigma_lmps_0, 3))))
	# ljcut.setParameter('sig',oxygenHandle, 'lmps_0', pow((pow(sigma, 6) + pow(sigma_lmps_0, 6))/2, 1/6))

 # 	ljcut.setParameter('sig',oxygenHandle, 'lmps_1', 2*math.sqrt(epsilon*epsilon_lmps_1)*((pow(sigma, 3) + pow(sigma_lmps_1, 3))/(pow(sigma, 6) + pow(sigma_lmps_1, 3))))
 # 	ljcut.setParameter('sig',oxygenHandle, 'lmps_1', pow((pow(sigma, 6) + pow(sigma_lmps_1, 6))/2, 1/6))

	# ljcut.setParameter('eps',oxygenHandle, 'lmps_2', 2*math.sqrt(epsilon*epsilon_lmps_2)*((pow(sigma, 3) + pow(sigma_lmps_2, 3))/(pow(sigma, 6) + pow(sigma_lmps_2, 3))))
	# ljcut.setParameter('sig',oxygenHandle, 'lmps_2', pow((pow(sigma, 6) + pow(sigma_lmps_2, 6))/2, 1/6))

	# ljcut.setParameter('eps',oxygenHandle, 'lmps_3', 2*math.sqrt(epsilon*epsilon_lmps_3)*((pow(sigma, 3) + pow(sigma_lmps_3, 3))/(pow(sigma, 6) + pow(sigma_lmps_3, 3))))
	# ljcut.setParameter('sig',oxygenHandle, 'lmps_3', pow((pow(sigma, 6) + pow(sigma_lmps_3, 6))/2, 1/6))

	# ljcut.setParameter('eps',oxygenHandle, 'lmps_4', 2*math.sqrt(epsilon*epsilon_lmps_4)*((pow(sigma, 3) + pow(sigma_lmps_4, 3))/(pow(sigma, 6) + pow(sigma_lmps_4, 3))))
	# ljcut.setParameter('sig',oxygenHandle, 'lmps_4', pow((pow(sigma, 6) + pow(sigma_lmps_4, 6))/2, 1/6))



	# Discover the bounds of the system in order to construct a tight box


	# Set an initial guess
	bounds_xlo = hexane_xlo
	bounds_xhi = hexane_xhi
	bounds_ylo = hexane_ylo
	bounds_yhi = hexane_yhi
	bounds_zlo = hexane_zlo
	bounds_zhi = hexane_zhi


	# Loop over the new unwrapped hexane coordinates
	for i in range(int(num_atoms_hexane[0])):

		x_pos = coordinates_hexane[i, 0]
		y_pos = coordinates_hexane[i, 1]
		z_pos = coordinates_hexane[i, 2]

		if x_pos > bounds_xhi:
			bounds_xhi = x_pos
		if x_pos < bounds_xlo:
			bounds_xlo = x_pos

		if y_pos > bounds_yhi:
			bounds_yhi = y_pos
		if y_pos < bounds_ylo:
			bounds_ylo = y_pos

		if z_pos > bounds_zhi:
			bounds_zhi = z_pos
		if z_pos < bounds_zlo:
			bounds_zlo = z_pos


	# Loop over the new translated and unwrapped water coordinates
	for i in range(3650):
		for j in range(4):

			x_pos = positions[i, j, 0]
			y_pos = positions[i, j, 1]
			z_pos = positions[i, j, 2]

			if x_pos > bounds_xhi:
				bounds_xhi = x_pos
			if x_pos < bounds_xlo:
				bounds_xlo = x_pos

			if y_pos > bounds_yhi:
				bounds_yhi = y_pos
			if y_pos < bounds_ylo:
				bounds_ylo = y_pos

			if z_pos > bounds_zhi:
				bounds_zhi = z_pos
			if z_pos < bounds_zlo:
				bounds_zlo = z_pos



	loVector = Vector(bounds_xlo, bounds_ylo, bounds_zlo)
	hiVector = Vector(bounds_xhi, bounds_yhi, bounds_zhi)


	state.bounds = Bounds(state, lo = loVector, hi = hiVector)




	dataTemp = state.dataManager.recordTemperature(interval=100)
	InitializeAtoms.initTemp(state, 'all', the_temp)

	dataBounds = state.dataManager.recordBounds(interval=100)



	# state.nPerRingPoly = nBeads
	# state.preparePIMD(the_temp)




	writeconfig = WriteConfig(state, handle='writer', fn='output_interface_classical_geometric', writeEvery=printFreq, format='xyz')
	state.activateWriteConfig(writeconfig)

	# writeRestart = WriteConfig(state, handle='restart', fn='output_interface_1bead*', format='xml', writeEvery=printFreq)
	# state.activateWriteConfig(writeRestart)

	writeconfig.write()

	#exit()







	integVerlet = IntegratorVerlet(state)
	integVerlet.run(nSteps)





def main(argv):

    parser = create_parser()
    args = parser.parse_args()
    files, options = convert_args(args)

    process_datafile(files)

if __name__ == "__main__":
    main(sys.argv[1:])


