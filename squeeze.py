from simtk.openmm.app import PDBFile, ForceField
from simtk.unit import *
import mdtools
import glob
#import matplotlib.pyplot as plt
#import seaborn as sns
#sns.set_context("notebook", font_scale=1.3)
import mdtraj

forcefield = ForceField('amber14/protein.ff14SB.xml', 
                        'amber14/tip3p.xml')
# Load and solvate PDB
pdb = PDBFile("pdz3_rat_apo_refined43_final.pdb")
mdsystem = mdtools.LatticeMDSystem(pdb.topology, pdb.positions, forcefield, "P 41 3 2")
mdsystem.buildSuperCell(1, 1, 1)
mdsystem.addSolvent(neutralize=True, positiveIon="Na+", negativeIon="Cl-")
mdsystem.save("prepped.pdb") # ! manually check if pdb file is correct
mdsystem.calmdown(posre=True)

# Squeeze -- 0.05% tolerance of target box volume. 
mdsystem.squeeze(tolerance=0.0005, maxIterations=20, maxSimtime=20*nanoseconds)
mdsystem.buildSimulation(filePrefix=f"squeezed_production1", 
                         saveTrajectory=True, trajInterval=50000, 
                         saveStateData=True, stateDataInterval=50000, ensemble="NVT")
mdsystem.saveCheckpoint("checkpoint"+str(0))
# intermediate checkpoints prevents losing stuff
for i in range(5):
    mdsystem.simulate(20*nanoseconds)
    mdsystem.saveCheckpoint("checkpoint"+str(i+1))


# Plot squeeze run
#h5trajs = sorted(glob.glob("iter*.h5"))

#sns.set_palette(sns.cubehelix_palette(len(h5trajs), start=.5, rot=-.75))
#plt.figure(figsize=(9, 6))
#for h5traj in h5trajs:
#    traj = mdtraj.load(h5traj)
#    n = int(len(traj.topology.select("water")) / 3)
#    plt.plot(traj.time, traj.unitcell_volumes, label=f"{n} waters")
#
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.xlabel("Time (ps)")
#plt.ylabel(r"Box Volume (nm$^3$)")
#plt.xlim(0, traj.time[-1])
#plt.tight_layout()
#plt.savefig("DHFR_321_squeeze2.png")

