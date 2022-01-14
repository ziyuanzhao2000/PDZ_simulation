# This is used to simulate a single PDZ domain 

from simtk.openmm.app import ForceField
from simtk.openmm.app import PDBFile
from simtk.unit import *

from mdtools import SolvatedMDSystem
import os, sys, getopt
import tqdm

# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "Dvsed:i:o:n:t:w:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)


forcefield = ForceField('amber/protein.ff14SB.xml', 
                        'amber14/tip3p.xml')

verbose = False
dryrun = False
inputname = "pdz3_rat_apo_refined43_final"
outputname = "production"
n_phases = 1000
t_per_phase = 10
t_eq = 10
box_width = 65 # angstrom
solvate = True
equilibriate = True

for o, a in opts:
    if o == "-d":
        dirname = a
    elif o == "-v":
        verbose = True
    elif o == "-D":
        dryrun = True
    elif o == "-i":
        inputname = a
    elif o == "-o":
        outputname = a
    elif o == "-n":
        n_phases = int(a)
    elif o == "-t":
        t_per_phase = int(a)
    elif o == "-s":
        solvate = False
    elif o == "-e":
        equilibriate = False
    elif o == "-w":
        box_width = float(a)

# Load PDB and create system
pdb = PDBFile(f"{inputname}.pdb")

if dirname != "":
    try:
        os.mkdir(dirname)
    except:
        pass # ignore err due to existing dir
    os.chdir(dirname)

mdsystem = SolvatedMDSystem(pdb.topology, pdb.positions, forcefield)

if solvate:
    mdsystem.addSolvent(boxSize=(box_width, box_width, box_width)*angstroms, positiveIon="Na+", negativeIon="Cl-")
    mdsystem.save("prepped.pdb")

if dryrun:
    sys.exit(1)

if equilibriate:
    mdsystem.buildSimulation(posre=True)
    mdsystem.equilibrate(t_eq*nanoseconds, posre=True)
    mdsystem.save("equilibrated.pdb")

# Production run
mdsystem.buildSimulation(filePrefix=f"{outputname}",
                         saveTrajectory=True, trajInterval=50000, 
                         saveStateData=True, stateDataInterval=50000)
mdsystem.saveCheckpoint("checkpoint_0")

for i in tqdm.tqdm(range(n_phases)):
    mdsystem.simulate(t_per_phase*nanoseconds)
    mdsystem.saveCheckpoint(f"checkpoint_{i+1}")
    if verbose:
        tqdm.tqdm.write(f"Phase {i} completed, simulated for {i * t_per_phase} ns in total!")



