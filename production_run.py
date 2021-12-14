from simtk.openmm.app import PDBFile, ForceField
from simtk.unit import *
import mdtools

import os
import sys, getopt

# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "Dvd:t:N:i:o:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

dirname = ""
maxit = 5
dt = 20
verbose = False
dryrun = False
inputname = "squeezed.pdb"
outputname = "production"

for o, a in opts:
    if o == "-d":
        dirname = a
    elif o == "-t":
        dt = int(a)
    elif o == "-N":
        maxit = int(a)
    elif o == "-v":
        verbose = True
    elif o == "-D":
        dryrun = True
    elif o == "-i":
        inputname = a
    elif o == "-o":
        outputname = a

if dirname != "":
    try: 
        os.mkdir(dirname)
    except:
        pass # ignore err due to existing dir 
    os.chdir(dirname)

pdb = PDBFile(inputname)

# Starts production run on the squeezed system
mdsystem.buildSimulation(filePrefix=outputname, 
                         saveTrajectory=True, trajInterval=50000,                       # 50000 timestep * 0.002 ps = one save every 100 ns
                         saveStateData=True, stateDataInterval=50000, ensemble="NVT")

if dryrun:
    sys.exit(1)

mdsystem.saveCheckpoint("checkpoint_"+str(0))

# Intermediate checkpoints prevents losing stuff
for iter in range(maxit):
    mdsystem.simulate(dt * nanoseconds)
    mdsystem.saveCheckpoint("checkpoint_" + str(iter+1))