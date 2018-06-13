"Script to run a DASH simulation of TIP4P/Flexible water, by Kirk Swanson"

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
sys.path = sys.path + ['/home/swansonk1/NEW_DASH_7.5/md_engine/build/python/build/lib.linux-x86_64-2.7']
sys.path.append('/home/swansonk1/NEW_DASH_7.5/md_engine/util_py')

# Import DASH packages
from DASH import *
from LAMMPS_Reader import LAMMPS_Reader
import re
import water
from water import *

# Functions to parse arguments
def create_parser():

    # Create parser and add arguments
    parser = argparse.ArgumentParser(description='Read data files')
    parser.add_argument('-equil_ensemble', dest='equil_ensemble', default='NPT', help='Thermodynamic ensemble for equilibration steps')
    parser.add_argument('-nSteps_equilibration', dest='nSteps_equilibration', default=0, help='Number of equilibration steps')
    parser.add_argument('-prod_ensemble', dest='prod_ensemble', default='NVT', help='Thermodynamic ensemble for production steps')
    parser.add_argument('-nSteps_production', dest='nSteps_production', default=0, help='Number of production steps')
    parser.add_argument('-record_restart', dest='record_restart', default=False, help='Boolean for recording xml restart files')
    parser.add_argument('-restartFreq', dest='restartFreq', default=10000, help='Frequency at which restart files are recorded')
    parser.add_argument('-record_traj', dest='record_traj', default=True, help='Boolean for recording configurations')
    parser.add_argument('-trajFreq1', dest='trajFreq1', default=1000, help='One frequency at which configurations are recorded')
    parser.add_argument('-trajFreq2', dest='trajFreq2', default=10000, help='A second frequency at which configurations are recorded')
    parser.add_argument('-numMolecules', dest='numMolecules', default=3650, help='Number of water molecules')
    parser.add_argument('-T', dest='T', default=298.15, help='Simulation temperature')
    parser.add_argument('-P', dest='P', default=1.0, help='Simulation pressure')
    parser.add_argument('-PI', dest='PI', default=True, help='Boolean for performing a path integral simulation')
    parser.add_argument('-nBeads', dest='nBeads', default=1, help='Number of beads for path integral simulation')
    parser.add_argument('-dataFreq', dest='dataFreq', default=20, help='Frequency at which pressure and pressure tensors are recorded')
    parser.add_argument('-zlo_change', dest='zlo_change', default=0.0, help='Change to the zlo bound')
    parser.add_argument('-zhi_change', dest='zhi_change', default=0.0, help='Change to the zhi bound')
    parser.add_argument('-filename', dest='filename', default='output_interface', help='Name of output files')
    parser.add_argument('-restart', dest='restart', default=False, help='Boolean for restarting from a previous xml configuration file')
    parser.add_argument('-restart_file', dest='restart_file', default=None, help='Name of DASH restart xml file')

    return parser

def convert_args(args):

    # Files dictionary
    files={}
    files['restart_file'] = args.restart_file

    options = {}
    options['equil_ensemble'] = args.equil_ensemble
    options['nSteps_equilibration'] = args.nSteps_equilibration
    options['prod_ensemble'] = args.prod_ensemble
    options['nSteps_production'] = args.nSteps_production
    options['record_restart'] = args.record_restart
    options['restartFreq'] = args.restartFreq
    options['record_traj'] = args.record_traj
    options['trajFreq1'] = args.trajFreq1
    options['trajFreq2'] = args.trajFreq2
    options['numMolecules'] = args.numMolecules
    options['T'] = args.T
    options['P'] = args.P
    options['PI'] = args.PI
    options['nBeads'] = args.nBeads
    options['dataFreq'] = args.dataFreq
    options['zlo_change'] = args.zlo_change
    options['zhi_change'] = args.zhi_change
    options['filename'] = args.filename
    options['restart'] = args.restart

    return files, options

""" Function to convert a string to a boolean for the argparse options"""
def str2bool(string):
    if string.lower() == 'true':
        return True
    elif string.lower() == 'false':
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# Main processing function
def process_datafile(files, options):

    # Parse options
    equil_ensemble = str(options['equil_ensemble'])
    nSteps_equilibration = int(options['nSteps_equilibration'])
    prod_ensemble = str(options['prod_ensemble'])
    nSteps_production = int(options['nSteps_production'])
    record_restart = str2bool(options['record_restart'])
    restartFreq = int(options['restartFreq'])
    record_traj = str2bool(options['record_traj'])
    trajFreq1 = int(options['trajFreq1'])
    trajFreq2 = int(options['trajFreq2'])
    numMolecules = int(options['numMolecules'])
    T = float(options['T'])
    P = float(options['P'])
    PI = str2bool(options['PI'])
    nBeads = int(options['nBeads'])
    dataFreq = int(options['dataFreq'])
    zlo_change = float(options['zlo_change'])
    zhi_change = float(options['zhi_change'])
    filename = str(options['filename'])
    restart = str2bool(options['restart'])

    if restart:

        # Initialize DASH settings
        state = State()
        state.deviceManager.setDevice(0)
        state.rCut = 9.0
        state.padding = 2.0
        state.periodicInterval = 1
        state.shoutEvery = 100
        state.units.setReal()
        state.dt = 0.500

        # Set the temperature
        the_temp = T

        # Adjust the temperature for path integral simulation if specified
        if PI:
            the_temp *= nBeads

        # Set atomic type labels
        oxygenHandle = 'OW'
        hydrogenHandle = 'HY'
        mSiteHandle = 'M'

        # Load restart file
        state.readConfig.loadFile(files['restart_file'])

        # Iterate to the first configuration
        state.readConfig.next()

        # LJ parameters
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
        fixNVT = FixNVTAndersen(state, handle='nvt', groupHandle='all', temp=the_temp, nu=0.01)
        fixNVT.setParameters(10)
        fixPressure = FixPressureBerendsen(state, 'npt', P, 1000, 5)

        # Activate fixes
        state.activateFix(nonbond)
        state.activateFix(charge)
        state.activateFix(bondQuart)
        state.activateFix(flexibleTIP4P)
        state.activateFix(harmonicAngle)
        state.activateFix(fixNVT)
        state.activateFix(fixPressure)

        # Data recording
        tempData = state.dataManager.recordTemperature('all', interval = dataFreq)
        boundsData = state.dataManager.recordBounds(interval=dataFreq)
        dataPressure = state.dataManager.recordPressure(handle='all', mode='scalar', interval=dataFreq)
        dataPressureTensor = state.dataManager.recordPressure(handle='all', mode='tensor', interval=dataFreq)

        # Initialize the system at the specified temperature
        InitializeAtoms.initTemp(state, 'all', the_temp)

        # Initialize a path integral simulation if specified
        if PI:
            state.nPerRingPoly = nBeads

        # Record restart xml files if specified
        if record_restart:
            writeRestart = WriteConfig(state, handle='restart', fn=filename+str('*'), format='xml', writeEvery=restartFreq)
            state.activateWriteConfig(writeRestart)

        # Add the integrator
        integVerlet = IntegratorVerlet(state)

        # Warning and exit if path integrals are combined with NPT simulation:
        if PI and (equil_ensemble == 'NPT' or prod_ensemble == 'NPT'):

            print('Path integral simulations can only be run in the NVT ensemble')
            exit()

        # Compute the total system mass for density calculation
        mtot = 0.0
        for atom in state.atoms:
            mtot += atom.mass

        # Equilibration run
        if nSteps_equilibration > 0:

            # Record trajectory if specified
            if record_traj:
                writeconfig1 = WriteConfig(state, handle='writer1', fn=filename+'-equil1-freq'+str(trajFreq1), writeEvery=trajFreq1, format='xyz')
                state.activateWriteConfig(writeconfig1)
                # Write the first configuration
                writeconfig1.write()

                writeconfig2 = WriteConfig(state, handle='writer2', fn=filename+'-equil2-freq'+str(trajFreq2), writeEvery=trajFreq2, format='xyz')
                state.activateWriteConfig(writeconfig2)
                # Write the first configuration
                writeconfig2.write()

            # Set up appropriate ensemble and run
            if equil_ensemble == 'NVT':
                state.deactivateFix(fixPressure)
                integVerlet.run(nSteps_equilibration)
            elif equil_ensemble == 'NPT':
                state.activateFix(fixPressure)
                integVerlet.run(nSteps_equilibration)

            # Deactivate the equilibration writer
            if record_traj:
                state.deactivateWriteConfig(writeconfig1)
                state.deactivateWriteConfig(writeconfig2)

        # Production run
        if nSteps_production > 0:

            # Record trajectory if specified
            if record_traj:
                writeconfig3 = WriteConfig(state, handle='writer3', fn=filename+'-prod1-freq'+str(trajFreq1), writeEvery=trajFreq1, format='xyz')
                state.activateWriteConfig(writeconfig3)
                # Write the first configuration
                writeconfig3.write()

                writeconfig4 = WriteConfig(state, handle='writer4', fn=filename+'-prod2-freq'+str(trajFreq2), writeEvery=trajFreq2, format='xyz')
                state.activateWriteConfig(writeconfig4)
                # Write the first configuration
                writeconfig4.write()

            # Set up appropriate ensemble and run
            if prod_ensemble == 'NVT':
                state.deactivateFix(fixPressure)
                integVerlet.run(nSteps_production)
            elif prod_ensemble == 'NPT':
                state.activateFix(fixPressure)
                integVerlet.run(nSteps_production)

        # Print list of total pressures
        output = open('pressures-'+filename+'.txt','w')
        for pressure in dataPressure.vals:
            output.write(str(pressure) + '\n')
        output.close()

        # Print list of pressure tensor elements
        output = open('pressure_tensors-'+filename+'.txt','w')
        for tensors in dataPressureTensor.vals:
            output.write(str(tensors[0]) + '\t')
            output.write(str(tensors[1]) + '\t')
            output.write(str(tensors[2]) + '\n')
        output.close()

        # Print list of densities
        output = open('densities-'+filename+'.txt', 'w')
        for bounds in boundsData.vals:
            output.write(str((1.0/0.6022)*(mtot/bounds.volume())) + '\n')
        output.close()

        exit()


    # Initialize DASH settings
    state = State()
    state.deviceManager.setDevice(0)
    Mw           = 1.0074*2 + 15.9994
    density      = 0.997
    Na           = 6.022140857e23
    sideLength   = (numMolecules * Mw * 1e24 / density / Na )**(1.0/3.0)
    print("The side length is "+str(sideLength))
    loVector = Vector(-0.5*sideLength,-0.5*sideLength,-0.5*sideLength)
    hiVector = Vector(0.5*sideLength, 0.5*sideLength, 0.5*sideLength)
    state.bounds = Bounds(state, lo = loVector, hi = hiVector)
    state.units.setReal()
    state.rCut = 9.0
    state.padding = 2.0
    state.periodicInterval = 1
    state.shoutEvery = 100
    state.dt = 0.500

    # Set the temperature
    the_temp = T

    # Adjust the temperature for path integral simulation if specified
    if PI:
        the_temp *= nBeads

    # Add molecular species
    oxygenHandle = 'OW'
    hydrogenHandle = 'HY'
    mSiteHandle = 'M'
    state.atomParams.addSpecies(handle=oxygenHandle, mass=15.9994, atomicNum=8)
    state.atomParams.addSpecies(handle=hydrogenHandle, mass=1.0074, atomicNum=1)
    state.atomParams.addSpecies(handle=mSiteHandle,mass=0.0,atomicNum=0)

    # LJ parameters
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
    fixNVT = FixNVTAndersen(state, handle='nvt', groupHandle='all', temp=the_temp, nu=0.01)
    fixNVT.setParameters(10)
    fixPressure = FixPressureBerendsen(state, 'npt', P, 1000, 5)

    # Add molecules
    positions = []

    xyzrange = int(math.ceil(numMolecules**(1.0/3.0)))
    xyzFloat = float(xyzrange)
    for x in xrange(xyzrange):
        for y in xrange(xyzrange):
            for z in xrange(xyzrange):
                pos = Vector( float(x)/(xyzFloat)*(0.95*sideLength) - 0.475*sideLength,
                              float(y)/(xyzFloat)*(0.95*sideLength) - 0.475*sideLength,
                              float(z)/(xyzFloat)*(0.95*sideLength) - 0.475*sideLength)

                positions.append(pos)

    for i in range(numMolecules):
        center = positions[i]

        molecule = create_TIP4P_Flexible(state,oxygenHandle,hydrogenHandle,mSiteHandle,center,"random")

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

    # Activate fixes
    state.activateFix(nonbond)
    state.activateFix(charge)
    state.activateFix(bondQuart)
    state.activateFix(flexibleTIP4P)
    state.activateFix(harmonicAngle)
    state.activateFix(fixNVT)
    state.activateFix(fixPressure)

    # Data recording
    tempData = state.dataManager.recordTemperature('all', interval = dataFreq)
    boundsData = state.dataManager.recordBounds(interval=dataFreq)
    dataPressure = state.dataManager.recordPressure(handle='all', mode='scalar', interval=dataFreq)
    dataPressureTensor = state.dataManager.recordPressure(handle='all', mode='tensor', interval=dataFreq)

    # Adjust the initial bounds if specified
    loVector = Vector(-0.5*sideLength,-0.5*sideLength,-0.5*sideLength + zlo_change)
    hiVector = Vector(0.5*sideLength, 0.5*sideLength, 0.5*sideLength + zhi_change)
    state.bounds = Bounds(state, lo = loVector, hi = hiVector)

    # Initialize the system at the specified temperature
    InitializeAtoms.initTemp(state, 'all', the_temp)

    # Initialize a path integral simulation if specified
    if PI:
        state.nPerRingPoly = nBeads
        state.preparePIMD(the_temp)

    # Record restart xml files if specified
    if record_restart:
        writeRestart = WriteConfig(state, handle='restart', fn=filename+str('*'), format='xml', writeEvery=restartFreq)
        state.activateWriteConfig(writeRestart)


    # Add the integrator
    integVerlet = IntegratorVerlet(state)

    # Warning and exit if path integrals are combined with NPT simulation:
    if PI and (equil_ensemble == 'NPT' or prod_ensemble == 'NPT'):

        print('Path integral simulations can only be run in the NVT ensemble')
        exit()

    # Compute the total system mass for density calculation
    mtot = 0.0
    for atom in state.atoms:
        mtot += atom.mass

    # Equilibration run
    if nSteps_equilibration > 0:

        # Record trajectory if specified
        if record_traj:
            writeconfig1 = WriteConfig(state, handle='writer1', fn=filename+'-equil1-freq'+str(trajFreq1), writeEvery=trajFreq1, format='xyz')
            state.activateWriteConfig(writeconfig1)
            # Write the first configuration
            writeconfig1.write()

            writeconfig2 = WriteConfig(state, handle='writer2', fn=filename+'-equil2-freq'+str(trajFreq2), writeEvery=trajFreq2, format='xyz')
            state.activateWriteConfig(writeconfig2)
            # Write the first configuration
            writeconfig2.write()

        # Set up appropriate ensemble and run
        if equil_ensemble == 'NVT':
            state.deactivateFix(fixPressure)
            integVerlet.run(nSteps_equilibration)
        elif equil_ensemble == 'NPT':
            state.activateFix(fixPressure)
            integVerlet.run(nSteps_equilibration)

        # Deactivate the equilibration writer
        if record_traj:
            state.deactivateWriteConfig(writeconfig1)
            state.deactivateWriteConfig(writeconfig2)


    # Production run
    if nSteps_production > 0:

        # Record trajectory if specified
        if record_traj:
            writeconfig3 = WriteConfig(state, handle='writer3', fn=filename+'-prod1-freq'+str(trajFreq1), writeEvery=100, format='xyz')
            state.activateWriteConfig(writeconfig3)
            # Write the first configuration
            writeconfig3.write()

            writeconfig4 = WriteConfig(state, handle='writer4', fn=filename+'-prod2-freq'+str(trajFreq2), writeEvery=50, format='xyz')
            state.activateWriteConfig(writeconfig4)
            # Write the first configuration
            writeconfig4.write()

        # Set up appropriate ensemble and run
        if prod_ensemble == 'NVT':
            state.deactivateFix(fixPressure)
            integVerlet.run(nSteps_production)
        elif prod_ensemble == 'NPT':
            state.activateFix(fixPressure)
            integVerlet.run(nSteps_production)

    # Print list of total pressures
    output = open('pressures-'+filename+'.txt','w')
    for pressure in dataPressure.vals:
        output.write(str(pressure) + '\n')
    output.close()

    # Print list of pressure tensor elements
    output = open('pressure_tensors-'+filename+'.txt','w')
    for tensors in dataPressureTensor.vals:
        output.write(str(tensors[0]) + '\t')
        output.write(str(tensors[1]) + '\t')
        output.write(str(tensors[2]) + '\n')
    output.close()

    # Print list of densities
    output = open('densities-'+filename+'.txt', 'w')
    for bounds in boundsData.vals:
        output.write(str((1.0/0.6022)*(mtot/bounds.volume())) + '\n')
    output.close()

def main(argv):

    parser = create_parser()
    args = parser.parse_args()
    files, options = convert_args(args)

    process_datafile(files, options)

if __name__ == "__main__":
    main(sys.argv[1:])


