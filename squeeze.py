from simtk.openmm.app import PDBFile, ForceField
from simtk.unit import *
import mdtools
import os
import sys, getopt
# import mdtraj

# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "vd:t:i:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

dirname = ""
tol = 0.0005
maxit = 20
verbose = False

for o, a in opts:
    if o == "-d":
        dirname = a
    elif o == "-t":
        tol = float(a)
    elif o == "-i":
        maxit = int(a)
    elif o == "-v":
        verbose = True

forcefield = ForceField('amber14/protein.ff14SB.xml', 
                        'amber14/tip3p.xml')
# Load and solvate PDB
pdb = PDBFile("pdz3_rat_apo_refined43_final.pdb")

if dirname != "":
    os.mkdir(dirname) # ignore err due to existing dir
    os.chdir(dirname)

# Prepare the system
mdsystem = mdtools.LatticeMDSystem(pdb.topology, pdb.positions, forcefield, "P 41 3 2")
mdsystem.buildSuperCell(1, 1, 1)
mdsystem.addSolvent(neutralize=True, positiveIon="Na+", negativeIon="Cl-")
mdsystem.save("prepped.pdb") # Should manually check if the structure is correct before proceeding
if verbose:
    print('Unit super cell built, solvent added, system neutralized.')

mdsystem.calmdown(posre=True)
if verbose:
    print('System calmed down. Squeeze job about to begin.')

# Squeeze -- 0.05% tolerance of target box volume. 
mdsystem.squeeze(tolerance=tol, maxIterations=maxit, maxSimtime=maxit*nanoseconds)
mdsystem.save("squeezed.pdb")
if verbose:
    print('Squeeze job completed. MDSystem saved to squeezed.pdb file.')




