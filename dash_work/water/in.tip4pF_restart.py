"Script to restart a TIP4P/F water simulation from a DASH restart file in DASH, by Kirk Swanson"
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
import water
from water import *

# Function to parse arguments for one DASH restart file
def create_parser():

    # Create parser and add arguments
    parser = argparse.ArgumentParser(description='Read data files')
    parser.add_argument('-restart_file', dest='restart_file', default=None, help='Name of DASH restart xml file')
        
    return parser

# Function to convert arguments into a dictionary
def convert_args(args):

    # Files dictionary
    files={}
    files['restart_file'] = args.restart_file
        
    options = {}

    # Print confirmation
    print("**************************")
    print("# Data file: ", files['restart_file'])
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

    # Read restart file
    f1 = open(files['restart_file'], 'r')
    contents1 = f1.readlines()
    f1.close()

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

    # Calculate box dimensions
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
    
    
    # Initialize DASH
    state = State()
    state.deviceManager.setDevice(0)

    # Settings
    numMolecules = 3650
    nSteps       = 10000
    Mw           = 1.0074*2 + 15.9994
    density      = 0.997
    Na           = 6.022140857e23
    sideLength   = (numMolecules * Mw * 1e24 / density / Na )**(1.0/3.0)
    print("The side length is ", sideLength)

    loVector = Vector(-0.5*sideLength,-0.5*sideLength,-0.5*sideLength)
    hiVector = Vector(0.5*sideLength, 0.5*sideLength, 0.5*sideLength)

    state.units.setReal()

    state.bounds = Bounds(state, lo = loVector, hi = hiVector)
    state.rCut = 9.0
    state.padding = 2.0
    state.periodicInterval = 1
    state.shoutEvery = 100
    state.dt = 0.500

    the_temp = 298.15

    # Add molecular species
    oxygenHandle = 'OW'
    hydrogenHandle = 'HY'
    mSiteHandle = 'M'

    state.atomParams.addSpecies(handle=oxygenHandle, mass=15.9994, atomicNum=8)
    state.atomParams.addSpecies(handle=hydrogenHandle, mass=1.0074, atomicNum=1)
    state.atomParams.addSpecies(handle=mSiteHandle,mass=0.0,atomicNum=0)

    epsilon = 0.1852 # given in kcal/mol
    sigma = 3.1589 # given in Angstroms

    # Add fixes
    charge = FixChargeEwald(state, 'charge', 'all')
    charge.setError(0.01,state.rCut,3)


    nonbond = FixLJCut(state,'cut')
    nonbond.setParameter('sig',oxygenHandle, oxygenHandle, sigma)
    nonbond.setParameter('eps',oxygenHandle, oxygenHandle, epsilon)
    nonbond.setParameter('sig',hydrogenHandle, hydrogenHandle, 0.0)
    nonbond.setParameter('eps',hydrogenHandle, hydrogenHandle, 0.0)
    nonbond.setParameter('sig',mSiteHandle, mSiteHandle, 0.0)
    nonbond.setParameter('eps',mSiteHandle, mSiteHandle, 0.0)

    flexibleTIP4P = FixTIP4PFlexible(state,'TIP4PFlexible')

    bondQuart = FixBondQuartic(state,'bondQuart')
    bondQuart.setBondTypeCoefs(type=0,k2=607.19,k3=-1388.65,k4=1852.58,r0=0.9419)

    harmonicAngle = FixAngleHarmonic(state,'angleH')
    harmonicAngle.setAngleTypeCoefs(type=0, k=87.85, theta0=( (107.4/180.0) * pi))

    # Change positions array to float
    positions = positions.astype(float)

    # Add molecules
    for i in range(numMolecules):
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

        harmonicAngle.createAngle(state.atoms[ids[1]], state.atoms[ids[0]], state.atoms[ids[2]],type=0)

    print('done adding molecules to simulation')
    #############################################################
    # Initialization of potentials
    #############################################################

    ############################################
    # Intermolecular interactions: LJ & Charge
    ############################################
    # -- we defined the LJ interactions above
    state.activateFix(nonbond)
    state.activateFix(charge)

    ################################################################
    # Intramolecular interactions:
    # FixTIP4PFlexible, quartic bonds, and harmonic angle potential
    ################################################################
    state.activateFix(bondQuart)
    state.activateFix(flexibleTIP4P)
    state.activateFix(harmonicAngle)

    integVerlet = IntegratorVerlet(state)

    ########################################
    # Data recording
    ########################################
    tempData = state.dataManager.recordTemperature('all', interval = 1)
    boundsData = state.dataManager.recordBounds(interval=100)
    ################################################
    # Thermostatting
    ################################################
    fixNVT = FixNVTAndersen(state, handle='nvt', groupHandle='all', temp=the_temp, nu=0.01)
    fixNVT.setParameters(10)
    state.activateFix(fixNVT)
    fixPressure = FixPressureBerendsen(state, 'npt', 1.0, 10000, 5)
    state.activateFix(fixPressure)
    InitializeAtoms.initTemp(state, 'all', the_temp)


    # Write configurations
    writer = WriteConfig(state, handle='writer', fn='tip4pF_new', format='xyz',
                         writeEvery=1000)
    state.activateWriteConfig(writer)

    writer.write();

    # Compute mass
    mtot = 0.0
    for atom in state.atoms:
        mtot += atom.mass

    # Run simulation
    integVerlet.run(nSteps)

    # Compute densities
    densities = []
    for bounds in boundsData.vals:
        densities.append((1.0/0.6022)*(mtot/bounds.volume()))

    print(densities)

        
def main(argv):

    parser = create_parser()
    args = parser.parse_args()
    files, options = convert_args(args)

    process_datafile(files)

if __name__ == "__main__":
    main(sys.argv[1:])




