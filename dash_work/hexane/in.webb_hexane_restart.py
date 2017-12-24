import sys
import numpy as np

sys.path = sys.path + ['/home/swansonk1/OLD_DASH/md_engine/build/python/build/lib.linux-x86_64-2.7']
sys.path.append('/home/swansonk1/OLD_DASH/md_engine/util_py')

from DASH import *
from LAMMPS_Reader import LAMMPS_Reader
import argparse
import re
import matplotlib.pyplot as plt
from math import *

state = State()
state.deviceManager.setDevice(0)
#state.rCut = 10.0
#state.padding = 2.0
#state.periodicInterval = 7
state.shoutEvery = 100
state.units.setReal()
#state.dt = 1.0
nSteps   = 500000
the_temp = 298.15
printFreq= 1000


state.readConfig.loadFile('output_webb_hexane_preequil950000.xml')
state.readConfig.next()

ljcut = FixLJCut(state, 'ljcut')
bondHARM = FixBondHarmonic(state, 'bondHARM')
angleHARM  = FixAngleHarmonic(state, 'angleHARM')
dihedralOPLS = FixDihedralOPLS(state, 'dihedralOPLS')

state.activateFix(ljcut)
state.activateFix(bondHARM)
state.activateFix(angleHARM)
state.activateFix(dihedralOPLS)

ewald = FixChargeEwald(state, "ewald", "all")
ewald.setParameters(32, state.rCut, 3)
state.activateFix(ewald)

fixNVT = FixNVTAndersen(state, handle='nvt', groupHandle='all', temp=the_temp, nu=0.01)
fixNVT.setParameters(10)
state.activateFix(fixNVT)

fixPressure = FixPressureBerendsen(state, 'npt', 1.0, 1000, 5)
state.activateFix(fixPressure)


writeconfig = WriteConfig(state, handle='writer', fn='output_webb_hexane_preequil_check', writeEvery=printFreq, format='xyz')
state.activateWriteConfig(writeconfig)

writeRestart = WriteConfig(state, handle='restart', fn='output_webb_hexane_preequil_check*', format='xml', writeEvery=50000)
state.activateWriteConfig(writeRestart)

#reader = LAMMPS_Reader(state=state, setBounds=True, nonbondFix = ljcut, bondFix = bondHARM, angleFix = angleHARM, dihedralFix = dihedralOPLS, atomTypePrefix = 'lmps_')
#reader.read(dataFn = 'in_webb_hexane.txt', inputFns = ['hexane.in.settings'])

dataTemp = state.dataManager.recordTemperature(interval=100)
InitializeAtoms.initTemp(state, 'all', the_temp)

dataBounds = state.dataManager.recordBounds(interval=100)

mtot = 0.0
for atom in state.atoms:
    mtot += atom.mass


integVerlet = IntegratorVerlet(state)
integVerlet.run(nSteps)

densities = []
for bounds in dataBounds.vals:
    densities.append((1.0/0.6022)*(mtot/bounds.volume()))

print densities
print np.mean(densities[len(densities)-100:])


