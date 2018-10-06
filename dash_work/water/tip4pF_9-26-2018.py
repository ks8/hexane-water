"Script to run a DASH simulation of q-TIP4P/F water using DASH-9-26-2018, by Kirk Swanson"
# Load modules and set system path
import sys
import argparse
import os
import matplotlib.pyplot as plt

sys.path = sys.path + ['/home/swansonk1/DASH-9-26-2018/md_engine/build/python/build/lib.linux-x86_64-2.7']
sys.path.append('/home/swansonk1/DASH-9-26-2018/md_engine/util_py')

from DASH import *
import water
from water import *
import math
import numpy as np
import time

# Functions to parse arguments
def create_parser():

    # Create parser and add arguments
    parser = argparse.ArgumentParser(description='Read data files')
    parser.add_argument('-numMolecules', dest='numMolecules', default=216, help='Number of molecules')
    parser.add_argument('-nSteps', dest='nSteps', default=200000, help='Number of simulation steps')
    parser.add_argument('-ensemble', dest='ensemble', default='NPT', help='Simulation ensemble')
    parser.add_argument('-barostat', dest='barostat', default='Berendsen', help='Type of barostat to use, either Berendsen or MonteCarlo')
    parser.add_argument('-T', dest='T', default=298.0, help='Simulation temperature')
    parser.add_argument('-P', dest='P', default=1.0, help='Simulation pressure')
    parser.add_argument('-density_init', dest='density_init', default=0.997, help='Initial configuration density')
    parser.add_argument('-rCut', dest='rCut', default=9.0, help='Cutoff radius for short-range interactions') 
    parser.add_argument('-timeStep', dest='timeStep', default=0.5, help='MD time step')        
    parser.add_argument('-PI', dest='PI', default=False, help='Boolean for path integral simulation')
    parser.add_argument('-nBeads', dest='nBeads', default=1, help='Number of beads for path integral simulation')
    parser.add_argument('-record_traj', dest='record_traj', default=True, help='Boolean for recording configurations')
    parser.add_argument('-trajFreq', dest='trajFreq', default=1000, help='Frequency at which configurations are recorded')
    parser.add_argument('-record_restart', dest='record_restart', default=False, help='Boolean for recording xml restart files')
    parser.add_argument('-restartFreq', dest='restartFreq', default=1000, help='Frequency at which restart files are recorded')
    parser.add_argument('-filename', dest='filename', default='test', help='Output filename')
    parser.add_argument('-dataFreq', dest='dataFreq', default=100, help='Data recording frequency')
    parser.add_argument('-restart', dest='restart', default=False, help='Boolean for restarting from a DASH restart xml file')
    parser.add_argument('-restart_file', dest='restart_file', default=None, help='Name of DASH restart xml file')

    return parser

def convert_args(args):

    # Files dictionary
    files={}
    files['restart_file'] = args.restart_file

    options = {}
    options['numMolecules'] = args.numMolecules
    options['nSteps'] = args.nSteps
    options['ensemble'] = args.ensemble
    options['barostat'] = args.barostat 
    options['T'] = args.T
    options['P'] = args.P
    options['density_init'] = args.density_init
    options['rCut'] = args.rCut
    options['timeStep'] = args.timeStep
    options['PI'] = args.PI
    options['nBeads'] = args.nBeads
    options['record_traj'] = args.record_traj
    options['trajFreq'] = args.trajFreq
    options['record_restart'] = args.record_restart
    options['restartFreq'] = args.restartFreq
    options['filename'] = args.filename
    options['dataFreq'] = args.dataFreq
    options['restart'] = args.restart

    return files, options

# Function to convert a string to a boolean for the argparse options
def str2bool(string):
    if string.lower() == 'true':
        return True
    elif string.lower() == 'false':
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# Main processing function
def process_datafile(files, options):

    numMolecules = int(options['numMolecules'])
    nSteps = int(options['nSteps'])
    ensemble = str(options['ensemble'])
    barostat = str(options['barostat'])
    T = float(options['T'])
    P = float(options['P'])
    density_init = float(options['density_init'])
    rCut = float(options['rCut'])
    timeStep = float(options['timeStep'])
    PI = str2bool(options['PI'])
    nBeads = int(options['nBeads'])
    record_traj = str2bool(options['record_traj'])
    trajFreq = int(options['trajFreq'])
    record_restart = str2bool(options['record_restart'])
    restartFreq = int(options['restartFreq'])
    filename = str(options['filename'])
    dataFreq = int(options['dataFreq'])
    restart = str2bool(options['restart'])

    # Initialize DASH
    state = State()
    state.deviceManager.setDevice(0)

    # Set some simulation settings
    if not restart:
        Mw           = 1.0074*2 + 15.9994
        density      = density_init
        Na           = 6.022140857e23
        sideLength   = (numMolecules * Mw * 1e24 / density / Na )**(1.0/3.0)
        print "The side length is ", sideLength
        loVector = Vector(-0.5*sideLength,-0.5*sideLength,-0.5*sideLength)
        hiVector = Vector(0.5*sideLength, 0.5*sideLength, 0.5*sideLength)
        state.bounds = Bounds(state, lo = loVector, hi = hiVector)

    state.units.setReal()
    state.rCut = rCut
    state.padding = 2.0
    state.periodicInterval = 1
    state.shoutEvery = 100
    state.dt = timeStep
    the_temp = T

    # Adjust the temperature for path integral simulation if specified
    if PI:
        the_temp *= nBeads

    # Add molecular species
    oxygenHandle = 'OW'
    hydrogenHandle = 'HY'
    mSiteHandle = 'M'

    if not restart:
        state.atomParams.addSpecies(handle=oxygenHandle, mass=15.9994, atomicNum=8)
        state.atomParams.addSpecies(handle=hydrogenHandle, mass=1.0074, atomicNum=1)
        state.atomParams.addSpecies(handle=mSiteHandle,mass=0.0,atomicNum=0)

    if restart:
        # Load restart file
        state.readConfig.loadFile(files['restart_file'])

        # Iterate to the first configuration
        state.readConfig.next()

    # Set LJ parameters
    epsilon = 0.1852 # given in kcal/mol
    sigma = 3.1589 # given in Angstroms

    # Add some fixes
    charge = FixChargeEwald(state, 'charge', 'all')
    charge.setError(0.01, state.rCut, 3)

    nonbond = FixLJCut(state,'cut')
    nonbond.setParameter('sig', oxygenHandle, oxygenHandle, sigma)
    nonbond.setParameter('eps', oxygenHandle, oxygenHandle, epsilon)
    nonbond.setParameter('sig', hydrogenHandle, hydrogenHandle, 0.0)
    nonbond.setParameter('eps', hydrogenHandle, hydrogenHandle, 0.0)
    nonbond.setParameter('sig', mSiteHandle, mSiteHandle, 0.0)
    nonbond.setParameter('eps', mSiteHandle, mSiteHandle, 0.0)

    flexibleTIP4P = FixTIP4PFlexible(state, 'TIP4PFlexible')

    bondQuart = FixBondQuartic(state, 'bondQuart')
    bondQuart.setBondTypeCoefs(type=0, k2=607.19, k3=-1388.65, k4=1852.58, r0=0.9419)

    harmonicAngle = FixAngleHarmonic(state, 'angleH')
    harmonicAngle.setAngleTypeCoefs(type=0, k=87.85, theta0=( (107.4/180.0) * pi))

    rppot = FixRingPolyPot(state, "rppot", "all")

    if not restart:
        # Add the molecules
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

            molecule = create_TIP4P_Flexible(state, oxygenHandle, hydrogenHandle, mSiteHandle, center, "random")

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

        print 'done adding molecules to simulation'

    # Activate some fixes
    state.activateFix(nonbond)
    state.activateFix(charge)
    state.activateFix(bondQuart)
    state.activateFix(flexibleTIP4P)
    state.activateFix(harmonicAngle)
    state.activateFix(rppot)

    # Initialize integrator
    integVerlet = IntegratorVerlet(state)

    # Establish data recording
    tempData = state.dataManager.recordTemperature('all', interval = dataFreq)
    #dataPressureScalar = state.dataManager.recordPressure(handle='all', mode='scalar', interval=dataFreq)
    boundsData = state.dataManager.recordBounds(interval=dataFreq)
    engData = state.dataManager.recordEnergy(handle='all', mode='scalar', interval=dataFreq)
    #kinengData = state.dataManager.recordTemperature('all', mode='vector', interval = dataFreq)

    # Thermostats and barostats
    fixNVT = FixNVTAndersen(state, handle='nvt', groupHandle='all', temp=the_temp, nu=0.01)
    fixNVT.setParameters(10)
    state.activateFix(fixNVT)

    if ensemble == 'NPT':
        if barostat == 'Berendsen':
            print('Using Berendsen Barostat')
            fixPressure = FixPressureBerendsen(state, 'npt', P, 10000, 5)
        elif barostat == 'MonteCarlo':
            print('Using MonteCarlo Barostat')
            fixPressure = FixPressureMonteCarlo(state, handle='all', pressure=P, scale=0.002, applyEvery=25,tune=True,tuneFreq=200*25, mode='direct')
        
    # Initialize system temperature
    InitializeAtoms.initTemp(state, 'all', the_temp)

    # Initialize a path integral simulation if specified
    if PI:
        state.nPerRingPoly = nBeads

        if not restart:
            state.preparePIMD(the_temp)

    # Write configurations and restart files
    if record_traj:
        writer = WriteConfig(state, handle='writer', fn=filename, format='xyz', writeEvery=trajFreq)
        state.activateWriteConfig(writer)

    if record_restart:
        writeRestart = WriteConfig(state, handle = 'restart', fn=filename+str('_*'), format='xml',writeEvery=restartFreq)
        state.activateWriteConfig(writeRestart)

    # Compute mass
    mtot = 0.0
    for atom in state.atoms:
        mtot += atom.mass

    if PI:
        mtot /= nBeads

    # Create data files
    output = open('simulation_data-'+filename+'.txt', 'w')
    output.write(str('#')+'\"Potential Energy (kJ/mole)\",'+'\"Temperature (K)\",'+'\"Density (g/mL)\"'+'\n')
    output.close()

    # Printing pressure information
    def data_recording(currentTurn):

        # # Compute total kinetic energy
        # totalkineng = 0.0

        # for val in kinengData.vals[-1]:
        #     totalkineng += val 

        # Print data
        output = open('simulation_data-'+filename+'.txt', 'a')
        #output.write(str(dataPressureScalar.vals[-1]) + ',')
        output.write(str(engData.vals[-1]*4.184) + ',')
        #output.write(str(totalkineng*(1/(7.242972*(10**25)))*(1.6022*(10**23))) + ',')
        output.write(str(tempData.vals[-1]) + ',')
        output.write(str((1.0/0.6022)*(mtot/boundsData.vals[-1].volume())) + '\n')
        output.close()

    # Construct pressure python operation
    dataOperation = PythonOperation(handle = 'dataOp', operateEvery = dataFreq, operation = data_recording, synchronous=True)
    state.activatePythonOperation(dataOperation)

    # Activate pressure fix last
    if ensemble == 'NPT':
        state.activateFix(fixPressure)

    # Run the simulation
    integVerlet.run(nSteps)

    # Write final trajectory and restart files
    if record_traj:
        writer.write()

    if record_restart:
        writeRestart.write()


def main(argv):

    parser = create_parser()
    args = parser.parse_args()
    files, options = convert_args(args)

    process_datafile(files, options)

if __name__ == "__main__":
    main(sys.argv[1:])
