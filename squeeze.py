from simtk.openmm.app import PDBFile, ForceField
from simtk.unit import *
import mdtools
import os
import sys, getopt
# import mdtraj

# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "Dvd:t:N:i:o:w:n:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

dirname = ""
tol = 0.0005
maxit = 20
verbose = False
dryrun = False
inputname = "pdz3_rat_apo_refined43_final.pdb"
outputname = "squeezed.pdb"
initial_water_perturb = 0
t_eq = 10
dn = 1000

for o, a in opts:
    if o == "-d":
        dirname = a
    elif o == "-t":
        tol = float(a)
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
    elif o == "-w":
        initial_water_perturb = int(a)
    elif o == "-n":
        dn = int(a)

forcefield = ForceField('amber14/protein.ff14SB.xml', 
                        'amber14/tip3p.xml')
# Load and solvate PDB
pdb = PDBFile(f"{inputname}.pdb")

if dirname != "":
    try: 
        os.mkdir(dirname)
    except:
        pass # ignore err due to existing dir 
    os.chdir(dirname)

# Prepare the system
mdsystem = mdtools.LatticeMDSystem(pdb.topology, pdb.positions, forcefield, "P 41 3 2")
mdsystem.buildSuperCell(1, 1, 1)
mdsystem.addSolvent(neutralize=True, positiveIon="Na+", negativeIon="Cl-")

if dryrun: 
    sys.exit(1)

mdsystem.save("prepped.pdb") # Should manually check if the structure is correct before proceeding
if verbose:
    print('Unit super cell built, solvent added, system neutralized.')

mdsystem.calmdown(posre=True)
if verbose:
    print('System calmed down. Squeeze job about to begin.')

# Squeeze -- 0.05% tolerance of target box volume. 
mdsystem.squeeze(tolerance=tol, maxIterations=maxit, maxSimtime=5*nanoseconds,
                 initial_water_perturb=initial_water_perturb, dn = dn)
mdsystem.save(outputname)
if verbose:
    print(f'Squeeze job completed. MDSystem saved to {outputname}.pdb file.')

mdsystem.equilibrate(t_eq*nanoseconds, posre=True)
mdsystem.save("equilibrated.pdb")
if verbose:
    print('System equilibrated. MDSystem saved to equilibrated.pdb file.')



