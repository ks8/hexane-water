"Script to run a DASH simulation of TIP4P/Flexible water, by Kirk Swanson"
import sys
import os
import matplotlib.pyplot as plt

sys.path = sys.path + ['/home/swansonk1/OLD_DASH/md_engine/build/python/build/lib.linux-x86_64-2.7']
sys.path.append('/home/swansonk1/OLD_DASH/md_engine/util_py')

from DASH import *
import water
from water import *
import math
import numpy as np

# Initialize DASH
state = State()
state.deviceManager.setDevice(0)

# Settings
numMolecules = 3650
nSteps       = 2000000
Mw           = 1.0074*2 + 15.9994
density      = 0.997
Na           = 6.022140857e23
sideLength   = (numMolecules * Mw * 1e24 / density / Na )**(1.0/3.0)
print "The side length is ", sideLength

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

print 'done adding molecules to simulation'
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

# Write configurations and restart files
writer = WriteConfig(state, handle='writer', fn='kirk_configPIMD', format='xyz',
                     writeEvery=1000)
state.activateWriteConfig(writer)

writeRestart = WriteConfig(state, handle = 'restart',fn="kirk_tip4p_restart*", format='xml',writeEvery=50000)
state.activateWriteConfig(writeRestart)

writer.write();

# Compute mass
mtot = 0.0
for atom in state.atoms:
    mtot += atom.mass

# Run the simulation
integVerlet.run(nSteps)

# Compute densities
densities = []
for bounds in boundsData.vals:
    densities.append((1.0/0.6022)*(mtot/bounds.volume()))

print densities



