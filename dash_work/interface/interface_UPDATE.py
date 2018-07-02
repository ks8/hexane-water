"Script to perform a DASH hexane-water interface simulation, by Kirk Swanson"

# Import python packages
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

# Set the DASH build to be used
sys.path = sys.path + ['/home/swansonk1/NEW_DASH/md_engine/build/python/build/lib.linux-x86_64-2.7']
sys.path.append('/home/swansonk1/NEW_DASH/md_engine/util_py')

# Import DASH packages
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
    parser.add_argument('-hexane_restart_file', dest='hexane_restart_file', default=None, help='Name of hexane LAMMPS input file')
    parser.add_argument('-hexane_settings_file', dest='hexane_settings_file', default=None, help='Name of hexane LAMMPS settings file')
    parser.add_argument('-equil_ensemble', dest='equil_ensemble', default='NPT', help='Thermodynamic ensemble for equilibration steps')    
    parser.add_argument('-nSteps_equilibration', dest='nSteps_equilibration', default=0, help='Number of equilibration steps')
    parser.add_argument('-prod_ensemble', dest='prod_ensemble', default='NVT', help='Thermodynamic ensemble for production steps')
    parser.add_argument('-nSteps_production', dest='nSteps_production', default=0, help='Number of production steps')
    parser.add_argument('-mix', dest='mix', default='standard', help='Mixing rule, which is either standard or waldman')
    parser.add_argument('-record_restart', dest='record_restart', default=False, help='Boolean for recording xml restart files')
    parser.add_argument('-restartFreq', dest='restartFreq', default=10000, help='Frequency at which restart files are recorded')
    parser.add_argument('-record_traj', dest='record_traj', default=True, help='Boolean for recording configurations')
    parser.add_argument('-trajFreq', dest='trajFreq', default=1000, help='Frequency at which configurations are recorded')
    parser.add_argument('-T', dest='T', default=298.15, help='Simulation temperature')
    parser.add_argument('-PI', dest='PI', default=True, help='Boolean for performing a path integral simulation')
    parser.add_argument('-nBeads', dest='nBeads', default=1, help='Number of beads for path integral simulation')
    parser.add_argument('-ptensorFreq', dest='ptensorFreq', default=20, help='Frequency at which pressure tensors are recorded')
    parser.add_argument('-zlo_change', dest='zlo_change', default=0.0, help='Change to the zlo bound')
    parser.add_argument('-zhi_change', dest='zhi_change', default=0.0, help='Change to the zhi bound')
    parser.add_argument('-filename', dest='filename', default='output_interface', help='Name of output files')
    parser.add_argument('-restart', dest='restart', default=False, help='Boolean for restarting from a previous xml configuration file')
    parser.add_argument('-restart_file', dest='restart_file', default=None, help='Name of DASH restart xml file')

    return parser

# Auxiliary Functions
""" Function to parse args """
def convert_args(args):

    # Files dictionary
    files={}
    files['water_restart_file'] = args.water_restart_file
    files['hexane_restart_file'] = args.hexane_restart_file
    files['hexane_settings_file'] = args.hexane_settings_file
    files['restart_file'] = args.restart_file 

    options = {}
    options['equil_ensemble'] = args.equil_ensemble
    options['nSteps_equilibration'] = args.nSteps_equilibration
    options['prod_ensemble'] = args.prod_ensemble
    options['nSteps_production'] = args.nSteps_production
    options['mix'] = args.mix
    options['record_restart'] = args.record_restart
    options['restartFreq'] = args.restartFreq
    options['record_traj'] = args.record_traj
    options['trajFreq'] = args.trajFreq
    options['T'] = args.T
    options['PI'] = args.PI
    options['nBeads'] = args.nBeads
    options['ptensorFreq'] = args.ptensorFreq
    options['zlo_change'] = args.zlo_change
    options['zhi_change'] = args.zhi_change
    options['filename'] = args.filename
    options['restart'] = args.restart


    # Print confirmation
    print("**************************")
    print("# Data file: ", files['water_restart_file'])
    print("# Data file: ", files['hexane_restart_file'])
    print("**************************")
    print(" ")

    return files, options

""" Function to convert a string to a boolean for the argparse options"""
def str2bool(string):
	if string.lower() == 'true':
		return True
	elif string.lower() == 'false':
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

""" Function to identify indices containing a particular substring in a file """
def indices_containing_substring(the_list, substring):
    indices = []
    for i, s in enumerate(the_list):
        if substring in s:
            indices.append(i)
    return indices

# Main processing function
def process_datafile(files, options):

	# Parse options
	equil_ensemble = str(options['equil_ensemble'])
	nSteps_equilibration = int(options['nSteps_equilibration'])
	prod_ensemble = str(options['prod_ensemble'])
	nSteps_production = int(options['nSteps_production'])
	mix = str(options['mix'])
	record_restart = str2bool(options['record_restart'])
	restartFreq = int(options['restartFreq'])
	record_traj = str2bool(options['record_traj'])
	trajFreq = int(options['trajFreq'])
	T = float(options['T'])
	PI = str2bool(options['PI'])
	nBeads = int(options['nBeads'])
	ptensorFreq = int(options['ptensorFreq'])
	zlo_change = float(options['zlo_change'])
	zhi_change = float(options['zhi_change'])
	filename = str(options['filename'])
	restart = str2bool(options['restart'])


	if restart:

		# Set up initial DASH specifications  
		state = State()
		state.deviceManager.setDevice(0)
		state.rCut = 12.0
		state.padding = 2.0
		state.periodicInterval = 1
		state.shoutEvery = 100
		state.units.setReal()
		state.dt = 0.5
		the_temp = T

		# Set atomic type labels 
		oxygenHandle = 'OW'
		hydrogenHandle = 'HY'
		mSiteHandle = 'M'

		# Adjust the temperature for path integral simulation if specified 
		if PI:
			the_temp *= nBeads

		# Load restart file
		state.readConfig.loadFile(files['restart_file'])

		# Iterate to the last (most recent) configuration in the xml file
		state.readConfig.prev()

		# Add fixes 
		ljcut = FixLJCut(state, 'ljcut')
		bondHARM = FixBondHarmonic(state, 'bondHARM')
		angleHARM  = FixAngleHarmonic(state, 'angleHARM')
		dihedralOPLS = FixDihedralOPLS(state, 'dihedralOPLS')
		bondQuart = FixBondQuartic(state,'bondQuart')
		flexibleTIP4P = FixTIP4PFlexible(state,'TIP4PFlexible')
		charge = FixChargeEwald(state, 'charge', 'all')
		fixNVT = FixNVTAndersen(state, handle='nvt', groupHandle='all', temp=the_temp, nu=0.01)
		fixPressure = FixPressureBerendsen(state, 'npt', 1.0, 1000, 5)

		# Activate fixes
		state.activateFix(ljcut)	
		state.activateFix(bondHARM)
		state.activateFix(angleHARM)
		state.activateFix(dihedralOPLS)
		state.activateFix(bondQuart)
		state.activateFix(flexibleTIP4P)
		state.activateFix(charge)
		state.activateFix(fixNVT)
		state.activateFix(fixPressure)

		# Set fix specifications (epsilon in kcal/mol, sigma in Angstroms)
		water_epsilon = 0.1852 # given in kcal/mol
		water_sigma = 3.1589 # given in Angstroms

		ljcut.setParameter('sig',oxygenHandle, oxygenHandle, water_sigma)
		ljcut.setParameter('eps',oxygenHandle, oxygenHandle, water_epsilon)
		ljcut.setParameter('sig',hydrogenHandle, hydrogenHandle, 0.0)
		ljcut.setParameter('eps',hydrogenHandle, hydrogenHandle, 0.0)
		ljcut.setParameter('sig',mSiteHandle, mSiteHandle, 0.0)
		ljcut.setParameter('eps',mSiteHandle, mSiteHandle, 0.0)
		angleHARM.setAngleTypeCoefs(type=0, k=87.85, theta0=( (107.4/180.0) * pi))
		bondQuart.setBondTypeCoefs(type=0,k2=607.19,k3=-1388.65,k4=1852.58,r0=0.9419)
		charge.setError(0.01,state.rCut,3)
		fixNVT.setParameters(10)

		# Set LJ parameters for atom types from the hexane LAMMPS file 
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

		# Establish mixing rules, finishing ljcut fix specifications
		if mix == 'standard':


			ljcut.setParameter('eps',oxygenHandle, 'lmps_0', math.sqrt(water_epsilon*epsilon_lmps_0))
			ljcut.setParameter('sig',oxygenHandle, 'lmps_0', (water_sigma + sigma_lmps_0)/2)

			ljcut.setParameter('eps',oxygenHandle, 'lmps_1', math.sqrt(water_epsilon*epsilon_lmps_1))
			ljcut.setParameter('sig',oxygenHandle, 'lmps_1', (water_sigma + sigma_lmps_1)/2)

			ljcut.setParameter('eps',oxygenHandle, 'lmps_2', math.sqrt(water_epsilon*epsilon_lmps_2))
			ljcut.setParameter('sig',oxygenHandle, 'lmps_2', (water_sigma + sigma_lmps_2)/2)

			ljcut.setParameter('eps',oxygenHandle, 'lmps_3', math.sqrt(water_epsilon*epsilon_lmps_3))
			ljcut.setParameter('sig',oxygenHandle, 'lmps_3', (water_sigma + sigma_lmps_3)/2)

			ljcut.setParameter('eps',oxygenHandle, 'lmps_4', math.sqrt(water_epsilon*epsilon_lmps_4))
			ljcut.setParameter('sig',oxygenHandle, 'lmps_4', (water_sigma + sigma_lmps_4)/2)

		if mix == 'waldman':


			ljcut.setParameter('eps',oxygenHandle, 'lmps_0', 2*math.sqrt(water_epsilon*epsilon_lmps_0)*((pow(water_sigma, 3)*pow(sigma_lmps_0, 3))/(pow(water_sigma, 6) + pow(sigma_lmps_0, 6))))
			ljcut.setParameter('sig',oxygenHandle, 'lmps_0', pow((pow(water_sigma, 6) + pow(sigma_lmps_0, 6))/2, 1/6))

		 	ljcut.setParameter('eps',oxygenHandle, 'lmps_1', 2*math.sqrt(water_epsilon*epsilon_lmps_1)*((pow(water_sigma, 3)*pow(sigma_lmps_1, 3))/(pow(water_sigma, 6) + pow(sigma_lmps_1, 6))))
		 	ljcut.setParameter('sig',oxygenHandle, 'lmps_1', pow((pow(water_sigma, 6) + pow(sigma_lmps_1, 6))/2, 1/6))

			ljcut.setParameter('eps',oxygenHandle, 'lmps_2', 2*math.sqrt(water_epsilon*epsilon_lmps_2)*((pow(water_sigma, 3)*pow(sigma_lmps_2, 3))/(pow(water_sigma, 6) + pow(sigma_lmps_2, 6))))
			ljcut.setParameter('sig',oxygenHandle, 'lmps_2', pow((pow(water_sigma, 6) + pow(sigma_lmps_2, 6))/2, 1/6))

			ljcut.setParameter('eps',oxygenHandle, 'lmps_3', 2*math.sqrt(water_epsilon*epsilon_lmps_3)*((pow(water_sigma, 3)*pow(sigma_lmps_3, 3))/(pow(water_sigma, 6) + pow(sigma_lmps_3, 6))))
			ljcut.setParameter('sig',oxygenHandle, 'lmps_3', pow((pow(water_sigma, 6) + pow(sigma_lmps_3, 6))/2, 1/6))

			ljcut.setParameter('eps',oxygenHandle, 'lmps_4', 2*math.sqrt(water_epsilon*epsilon_lmps_4)*((pow(water_sigma, 3)*pow(sigma_lmps_4, 3))/(pow(water_sigma, 6) + pow(sigma_lmps_4, 6))))
			ljcut.setParameter('sig',oxygenHandle, 'lmps_4', pow((pow(water_sigma, 6) + pow(sigma_lmps_4, 6))/2, 1/6))


		# Initialize the system at the specified temperature
		InitializeAtoms.initTemp(state, 'all', the_temp)

		# Initialize a path integral simulation if specified 
		if PI:
			state.nPerRingPoly = nBeads

		# Record trajectory if specified
		if record_traj:
			writeconfig = WriteConfig(state, handle='writer', fn=filename, writeEvery=trajFreq, format='xyz')
			state.activateWriteConfig(writeconfig)
			

		# Record restart xml files if specified
		if record_restart:
			writeRestart = WriteConfig(state, handle='restart', fn=filename+str('_*'), format='xml', writeEvery=restartFreq)
			state.activateWriteConfig(writeRestart)

		# Add the integrator
		integVerlet = IntegratorVerlet(state)

		# Warning and exit if path integrals are combined with NPT simulation:
		if PI and (equil_ensemble == 'NPT' or prod_ensemble == 'NPT'):

			print('Path integral simulations can only be run in the NVT ensemble')
			exit()

		# Equilibration run 
		if nSteps_equilibration > 0:

			# Set up appropriate ensemble and run 
			if equil_ensemble == 'NVT':			
				state.deactivateFix(fixPressure)
				integVerlet.run(nSteps_equilibration)
			elif equil_ensemble == 'NPT':
				state.activateFix(fixPressure)
				integVerlet.run(nSteps_equilibration)

		# Production run 
		if nSteps_production > 0:

			# Record pressure
			dataPressureTensor = state.dataManager.recordPressure(handle='all', mode='tensor', interval=ptensorFreq)

			# Set up appropriate ensemble and run 
			if prod_ensemble == 'NVT':
				state.deactivateFix(fixPressure)
				integVerlet.run(nSteps_production)
			elif prod_ensemble == 'NPT':
				state.activateFix(fixPressure)
				integVerlet.run(nSteps_production)

			# Print list of pressures
			output = open('pressure_tensors-'+filename+'.txt','w')
			for tensors in dataPressureTensor.vals:
				output.write(str(tensors[0]) + '\t')
				output.write(str(tensors[1]) + '\t')
				output.write(str(tensors[2]) + '\n')
			output.close()

		# Write final configurations
		if record_traj:
			writeconfig.write()

		if record_restart:
			writeRestart.write()

		exit()






	# Read pure component restart files
	f1 = open(files['water_restart_file'], 'r')
	contents1 = f1.readlines()
	f1.close()

	f2 = open(files['hexane_restart_file'], 'r')
	contents_hexane = f2.readlines()
	f2.close()

	# Extract the simulation box coordinates from the hexane file
	hexane_x_bounds = contents_hexane[indices_containing_substring(contents_hexane, "xlo")[0]].split()
	hexane_xlo = float(hexane_x_bounds[0])
	hexane_xhi = float(hexane_x_bounds[1])

	hexane_y_bounds = contents_hexane[indices_containing_substring(contents_hexane, "ylo")[0]].split()
	hexane_ylo = float(hexane_y_bounds[0])
	hexane_yhi = float(hexane_y_bounds[1])

	hexane_z_bounds = contents_hexane[indices_containing_substring(contents_hexane, "zlo")[0]].split()
	hexane_zlo = float(hexane_z_bounds[0])
	hexane_zhi = float(hexane_z_bounds[1])

	# Save the LAMMPS header of the hexane file
	num_atoms_hexane = contents_hexane[2].split()
	num_bonds_hexane = contents_hexane[3].split()
	num_angles_hexane = contents_hexane[4].split()
	num_dihedrals_hexane = contents_hexane[5].split()

	num_atomtype_hexane = contents_hexane[7].split()
	num_bondtype_hexane = contents_hexane[8].split()
	num_angletype_hexane = contents_hexane[9].split()
	num_dihedraltype_hexane = contents_hexane[10].split()

	# Save the atomic coordinates of the hexane file
	coordinates_hexane = np.zeros([int(num_atoms_hexane[0]), 3])
	atom_index_hexane = contents_hexane.index("Atoms\n")
	for i in range(int(num_atoms_hexane[0])):
		coordinates_contents = contents_hexane[atom_index_hexane + 2 + i].split()
		coordinates_hexane[i][0] = float(coordinates_contents[4])
		coordinates_hexane[i][1] = float(coordinates_contents[5])
		coordinates_hexane[i][2] = float(coordinates_contents[6])

	# Set the number of atoms per molecule for hexane
	atoms_per_molecule_hexane = 20

	# Calculate the maximum dimensional range of the molecule (to prevent overlapping using appropriate spacing in between the initial hexane and water slabs along the z direction)
	ranges_hexane = max(np.max(coordinates_hexane[0:atoms_per_molecule_hexane], axis = 0) - np.min(coordinates_hexane[0:atoms_per_molecule_hexane], axis = 0))

	# Calculate number of atoms in the water restart file
	num_atoms_line = contents1[indices_containing_substring(contents1, "numAtoms")[0]].split()
	index = indices_containing_substring(contents1[indices_containing_substring(contents1, "numAtoms")[0]].split(), "numAtoms")[0]
	start = indices_containing_substring(num_atoms_line[index], '"')[0]
	end = indices_containing_substring(num_atoms_line[index], '"')[1]
	num_atoms = int(num_atoms_line[index][start+1:end])

	# Find the line where atomic positions are recorded in the water restart file
	position_line_index = indices_containing_substring(contents1, "<position>")[0]

	# Set the number of atoms per molecule for water
	atoms_per_molecule = 4

	# Find the line where simulation box position is recorded in the water restart file
	bounds_line = contents1[indices_containing_substring(contents1, "bounds")[0]].split()

	# Extract x lo position
	index_xlo = indices_containing_substring(bounds_line, "xlo")[0]
	start = indices_containing_substring(bounds_line[index_xlo], '"')[0]
	end = indices_containing_substring(bounds_line[index_xlo], '"')[1]
	water_xlo = float(bounds_line[index_xlo][start+1:end])

	# Extract x hi position
	index_xhi = indices_containing_substring(bounds_line, "xhi")[0]
	start = indices_containing_substring(bounds_line[index_xhi], '"')[0]
	end = indices_containing_substring(bounds_line[index_xhi], '"')[1]
	water_xhi = float(bounds_line[index_xhi][start+1:end])

	# Extract y lo position
	index_ylo = indices_containing_substring(bounds_line, "ylo")[0]
	start = indices_containing_substring(bounds_line[index_ylo], '"')[0]
	end = indices_containing_substring(bounds_line[index_ylo], '"')[1]
	water_ylo = float(bounds_line[index_ylo][start+1:end])

	# Extract y hi position
	index_yhi = indices_containing_substring(bounds_line, "yhi")[0]
	start = indices_containing_substring(bounds_line[index_yhi], '"')[0]
	end = indices_containing_substring(bounds_line[index_yhi], '"')[1]
	water_yhi = float(bounds_line[index_yhi][start+1:end])

	# Extract z lo position
	index_zlo = indices_containing_substring(bounds_line, "zlo")[0]
	start = indices_containing_substring(bounds_line[index_zlo], '"')[0]
	end = indices_containing_substring(bounds_line[index_zlo], '"')[1]
	water_zlo = float(bounds_line[index_zlo][start+1:end])

	# Extract z hi position
	index_zhi = indices_containing_substring(bounds_line, "zhi")[0]
	start = indices_containing_substring(bounds_line[index_zhi], '"')[0]
	end = indices_containing_substring(bounds_line[index_zhi], '"')[1]
	water_zhi = float(bounds_line[index_zhi][start+1:end])

	# Compute the box lengths
	x_boxlength = water_xhi - water_xlo
	y_boxlength = water_yhi - water_ylo
	z_boxlength = water_zhi - water_zlo

	# List to hold water atomic positions
	coordinates_water = []

	# Add atomic water positions to list 
	for i in range(int(num_atoms/atoms_per_molecule)):

		molecule = []
		molecule.append(contents1[position_line_index+1+atoms_per_molecule*i].split())
		molecule.append(contents1[position_line_index+2+atoms_per_molecule*i].split())
		molecule.append(contents1[position_line_index+3+atoms_per_molecule*i].split())
		molecule.append(contents1[position_line_index+4+atoms_per_molecule*i].split())

		coordinates_water.append(molecule)

	coordinates_water = np.array(coordinates_water)

	# Move the water slab so that it doesn't overlap with the hexane slab and so that there is appropriate spacing between slabs in the z direction
	for i in range(int(num_atoms/atoms_per_molecule)):
		for j in range(atoms_per_molecule):
			coordinates_water[i,j,0] = str(float(coordinates_water[i,j,0]) + hexane_xlo - water_xlo)
			coordinates_water[i,j,1] = str(float(coordinates_water[i,j,1]) + hexane_ylo - water_ylo)
			coordinates_water[i,j,2] = str(float(coordinates_water[i,j,2]) - (water_zhi - hexane_zlo + ranges_hexane))

	# Unwrap the water molecules
	for i in range(int(num_atoms/atoms_per_molecule)):

		molecule = list(coordinates_water[i, :, :])
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

			coordinates_water[i, j, :] = string_contents


	# Set up initial DASH specifications  
	state = State()
	state.deviceManager.setDevice(0)
	state.rCut = 12.0
	state.padding = 2.0
	state.periodicInterval = 1
	state.shoutEvery = 100
	state.units.setReal()
	state.dt = 0.5
	the_temp = T

	# Set atomic type labels and add atom types 
	oxygenHandle = 'OW'
	hydrogenHandle = 'HY'
	mSiteHandle = 'M'
	state.atomParams.addSpecies(handle=oxygenHandle, mass=15.9994, atomicNum=8)
	state.atomParams.addSpecies(handle=hydrogenHandle, mass=1.0074, atomicNum=1)
	state.atomParams.addSpecies(handle=mSiteHandle,mass=0.0,atomicNum=0)

	# Adjust the temperature for path integral simulation if specified 
	if PI:
		the_temp *= nBeads

	# Add fixes 
	ljcut = FixLJCut(state, 'ljcut')
	bondHARM = FixBondHarmonic(state, 'bondHARM')
	angleHARM  = FixAngleHarmonic(state, 'angleHARM')
	dihedralOPLS = FixDihedralOPLS(state, 'dihedralOPLS')
	bondQuart = FixBondQuartic(state,'bondQuart')
	flexibleTIP4P = FixTIP4PFlexible(state,'TIP4PFlexible')
	charge = FixChargeEwald(state, 'charge', 'all')
	fixNVT = FixNVTAndersen(state, handle='nvt', groupHandle='all', temp=the_temp, nu=0.01)
	fixPressure = FixPressureBerendsen(state, 'npt', 1.0, 1000, 5)

	# Activate fixes
	state.activateFix(ljcut)	
	state.activateFix(bondHARM)
	state.activateFix(angleHARM)
	state.activateFix(dihedralOPLS)
	state.activateFix(bondQuart)
	state.activateFix(flexibleTIP4P)
	state.activateFix(charge)
	state.activateFix(fixNVT)
	state.activateFix(fixPressure)

	# Set fix specifications (epsilon in kcal/mol, sigma in Angstroms)
	water_epsilon = 0.1852 # given in kcal/mol
	water_sigma = 3.1589 # given in Angstroms

	ljcut.setParameter('sig',oxygenHandle, oxygenHandle, water_sigma)
	ljcut.setParameter('eps',oxygenHandle, oxygenHandle, water_epsilon)
	ljcut.setParameter('sig',hydrogenHandle, hydrogenHandle, 0.0)
	ljcut.setParameter('eps',hydrogenHandle, hydrogenHandle, 0.0)
	ljcut.setParameter('sig',mSiteHandle, mSiteHandle, 0.0)
	ljcut.setParameter('eps',mSiteHandle, mSiteHandle, 0.0)
	angleHARM.setAngleTypeCoefs(type=0, k=87.85, theta0=( (107.4/180.0) * pi))
	bondQuart.setBondTypeCoefs(type=0,k2=607.19,k3=-1388.65,k4=1852.58,r0=0.9419)
	charge.setError(0.01,state.rCut,3)
	fixNVT.setParameters(10)

	# Add TIP4P/F water molecules to the simulation
	coordinates_water = coordinates_water.astype(float)
	for i in range(3650):

		center = coordinates_water[i, :, :]
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

	# Add hexane molecules to the simulation
	reader = LAMMPS_Reader(state=state, setBounds=True, nonbondFix = ljcut, bondFix = bondHARM, angleFix = angleHARM, dihedralFix = dihedralOPLS, atomTypePrefix = 'lmps_')
	reader.read(dataFn = files['hexane_restart_file'], inputFns = [files['hexane_settings_file']])

	# Set LJ parameters for atom types from the hexane LAMMPS file 
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

	# Establish mixing rules, finishing ljcut fix specifications
	if mix == 'standard':


		ljcut.setParameter('eps',oxygenHandle, 'lmps_0', math.sqrt(water_epsilon*epsilon_lmps_0))
		ljcut.setParameter('sig',oxygenHandle, 'lmps_0', (water_sigma + sigma_lmps_0)/2)

		ljcut.setParameter('eps',oxygenHandle, 'lmps_1', math.sqrt(water_epsilon*epsilon_lmps_1))
		ljcut.setParameter('sig',oxygenHandle, 'lmps_1', (water_sigma + sigma_lmps_1)/2)

		ljcut.setParameter('eps',oxygenHandle, 'lmps_2', math.sqrt(water_epsilon*epsilon_lmps_2))
		ljcut.setParameter('sig',oxygenHandle, 'lmps_2', (water_sigma + sigma_lmps_2)/2)

		ljcut.setParameter('eps',oxygenHandle, 'lmps_3', math.sqrt(water_epsilon*epsilon_lmps_3))
		ljcut.setParameter('sig',oxygenHandle, 'lmps_3', (water_sigma + sigma_lmps_3)/2)

		ljcut.setParameter('eps',oxygenHandle, 'lmps_4', math.sqrt(water_epsilon*epsilon_lmps_4))
		ljcut.setParameter('sig',oxygenHandle, 'lmps_4', (water_sigma + sigma_lmps_4)/2)

	if mix == 'waldman':


		ljcut.setParameter('eps',oxygenHandle, 'lmps_0', 2*math.sqrt(water_epsilon*epsilon_lmps_0)*((pow(water_sigma, 3)*pow(sigma_lmps_0, 3))/(pow(water_sigma, 6) + pow(sigma_lmps_0, 6))))
		ljcut.setParameter('sig',oxygenHandle, 'lmps_0', pow((pow(water_sigma, 6) + pow(sigma_lmps_0, 6))/2, 1/6))

	 	ljcut.setParameter('eps',oxygenHandle, 'lmps_1', 2*math.sqrt(water_epsilon*epsilon_lmps_1)*((pow(water_sigma, 3)*pow(sigma_lmps_1, 3))/(pow(water_sigma, 6) + pow(sigma_lmps_1, 6))))
	 	ljcut.setParameter('sig',oxygenHandle, 'lmps_1', pow((pow(water_sigma, 6) + pow(sigma_lmps_1, 6))/2, 1/6))

		ljcut.setParameter('eps',oxygenHandle, 'lmps_2', 2*math.sqrt(water_epsilon*epsilon_lmps_2)*((pow(water_sigma, 3)*pow(sigma_lmps_2, 3))/(pow(water_sigma, 6) + pow(sigma_lmps_2, 6))))
		ljcut.setParameter('sig',oxygenHandle, 'lmps_2', pow((pow(water_sigma, 6) + pow(sigma_lmps_2, 6))/2, 1/6))

		ljcut.setParameter('eps',oxygenHandle, 'lmps_3', 2*math.sqrt(water_epsilon*epsilon_lmps_3)*((pow(water_sigma, 3)*pow(sigma_lmps_3, 3))/(pow(water_sigma, 6) + pow(sigma_lmps_3, 6))))
		ljcut.setParameter('sig',oxygenHandle, 'lmps_3', pow((pow(water_sigma, 6) + pow(sigma_lmps_3, 6))/2, 1/6))

		ljcut.setParameter('eps',oxygenHandle, 'lmps_4', 2*math.sqrt(water_epsilon*epsilon_lmps_4)*((pow(water_sigma, 3)*pow(sigma_lmps_4, 3))/(pow(water_sigma, 6) + pow(sigma_lmps_4, 6))))
		ljcut.setParameter('sig',oxygenHandle, 'lmps_4', pow((pow(water_sigma, 6) + pow(sigma_lmps_4, 6))/2, 1/6))



	# Set the simulation bounds
	# Set an initial guess to be the bounds of the hexane slab
	bounds_xlo = hexane_xlo
	bounds_xhi = hexane_xhi
	bounds_ylo = hexane_ylo
	bounds_yhi = hexane_yhi
	bounds_zlo = hexane_zlo
	bounds_zhi = hexane_zhi

	# Loop over the new unwrapped hexane coordinates to readjust the initial guess
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

	# Loop over the new translated and unwrapped water coordinates to readjust again
	for i in range(3650):
		for j in range(4):

			x_pos = coordinates_water[i, j, 0]
			y_pos = coordinates_water[i, j, 1]
			z_pos = coordinates_water[i, j, 2]

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

	# Adjust the initial bounds if specified 
	bounds_zlo = bounds_zlo + zlo_change
	bounds_zhi = bounds_zhi + zhi_change 	

	# Set the simulation bounds
	loVector = Vector(bounds_xlo, bounds_ylo, bounds_zlo)
	hiVector = Vector(bounds_xhi, bounds_yhi, bounds_zhi)
	state.bounds = Bounds(state, lo = loVector, hi = hiVector)

	# Initialize the system at the specified temperature
	InitializeAtoms.initTemp(state, 'all', the_temp)

	# Initialize a path integral simulation if specified 
	if PI:
		state.nPerRingPoly = nBeads
		state.preparePIMD(the_temp)

	# Record trajectory if specified
	if record_traj:
		writeconfig = WriteConfig(state, handle='writer', fn=filename, writeEvery=trajFreq, format='xyz')
		state.activateWriteConfig(writeconfig)

	# Record restart xml files if specified
	if record_restart:
		writeRestart = WriteConfig(state, handle='restart', fn=filename+str('_*'), format='xml', writeEvery=restartFreq)
		state.activateWriteConfig(writeRestart)

	# Add the integrator
	integVerlet = IntegratorVerlet(state)

	# Warning and exit if path integrals are combined with NPT simulation:
	if PI and (equil_ensemble == 'NPT' or prod_ensemble == 'NPT'):

		print('Path integral simulations can only be run in the NVT ensemble')
		exit()

	# Equilibration run 
	if nSteps_equilibration > 0:

		# Set up appropriate ensemble and run 
		if equil_ensemble == 'NVT':			
			state.deactivateFix(fixPressure)
			integVerlet.run(nSteps_equilibration)
		elif equil_ensemble == 'NPT':
			state.activateFix(fixPressure)
			integVerlet.run(nSteps_equilibration)

	# Production run 
	if nSteps_production > 0:

		# Record pressure
		dataPressureTensor = state.dataManager.recordPressure(handle='all', mode='tensor', interval=ptensorFreq)

		# Set up appropriate ensemble and run 
		if prod_ensemble == 'NVT':
			state.deactivateFix(fixPressure)
			integVerlet.run(nSteps_production)
		elif prod_ensemble == 'NPT':
			state.activateFix(fixPressure)
			integVerlet.run(nSteps_production)

		# Print list of pressures
		output = open('pressure_tensors-'+filename+'.txt','w')
		for tensors in dataPressureTensor.vals:
			output.write(str(tensors[0]) + '\t')
			output.write(str(tensors[1]) + '\t')
			output.write(str(tensors[2]) + '\n')
		output.close()

	# Write final configurations
	if record_traj:
		writeconfig.write()

	if record_restart:
		writeRestart.write()

	exit()





def main(argv):

    parser = create_parser()
    args = parser.parse_args()
    files, options = convert_args(args)

    process_datafile(files, options)

if __name__ == "__main__":
    main(sys.argv[1:])


