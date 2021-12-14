# This is used to simulate a single PDZ domain 

from simtk.openmm.app import ForceField
from simtk.openmm.app import PDBFile
from simtk.unit import *

from mdtools import SolvatedMDSystem

forcefield = ForceField('amber/protein.ff14SB.xml', 
                        'tip3p.xml')

# Load and solvate PDB
pdb = PDBFile("pdz3_rat_apo_refined43_final.pdb")
mdsystem = SolvatedMDSystem(pdb.topology, pdb.positions, forcefield)
a = b = c = 90.013
mdsystem.addSolvent(boxSize=(a, b, c)*angstroms, ionicStrength=0.1*molar)
mdsystem.save("prepped.pdb")
mdsystem.saveCheckpoint("prepped_checkpoint"+0)

# Equilibrate
mdsystem.buildSimulation(posre=True)
mdsystem.equilibrate(10*nanoseconds, posre=True)
mdsystem.save("equilibrated.pdb")

# Production run
mdsystem.buildSimulation(filePrefix=f"production1", 
                         saveTrajectory=True, trajInterval=50000, 
                         saveStateData=True, stateDataInterval=50000)
for i in range(5):
    mdsystem.simulate(20*nanoseconds)
    mdsystem.saveCheckpoint("checkpoint"+(i+1))



