# Post-processing of MD production run trajectories of crystal system
# The script is based on chain_rmsd.ipynb

import os, sys, getopt
import mdtraj
import matplotlib.pyplot as plt
import numpy as np
from utils import smartWrapProtein

# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "svd:i:o:n:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

verbose = False
dir_name = ""
input_name = "production"
output_name = "production"
dim = 0
skip = False

for o, a in opts:
    if o == "-v":
        verbose = True # do we even need this??
    elif o == "-d":
        dir_name = a
    elif o == "-i":
        input_name = a
    elif o == "-o":
        output_name = a
    elif o == "-s":
        skip = True

if dir_name != "":
    try:
        os.mkdir(dir_name)
    except:
        pass # ignore err due to existing dir
    os.chdir(dir_name)


target_traj = mdtraj.load(f"{input_name}.h5") # currently only support h5 file!

if skip == False:
    smartWrapProtein(target_traj)
    target_traj.save_hdf5(f"{input_name}.h5") # overwrites
    if verbose:
        print("Wrapped proteins at box boundaries and saved the new trajectory.")

top = target_traj.topology

n_frames = target_traj.n_frames
unit_cell_rmsd = mdtraj.rmsd(target_traj, target_traj, 0, top.select("is_backbone"))
plt.rc('font', size=16)
fig1 = plt.figure(figsize=(18, 12))
ax = fig1.add_subplot()
ax.set_title("PDZ domain, unit cell RMSD")
ax.set_xlabel("Time (ns)")
ax.set_xticks(np.arange(0, n_frames,100))
ax.set_xticklabels(np.arange(0, n_frames,100) / 10)
ax.set_ylabel("RMSD (nm)")
ax.plot(unit_cell_rmsd, label = "unit cell 1")
ax.legend(loc = 'upper right')
fig1.save(f'{output_name}_unitcell_rmsd.pdf')
if verbose:
    print("Produced plot for unit-cell backbone RMSD from the trajectory.")