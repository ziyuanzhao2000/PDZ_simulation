# This script reads in the tICA-transformed trajectory stored in the ttrajs folder and then plot the time trace
# along each direction. But how to interpret this??

import numpy as np
from msmbuilder.io import load_trajs, load_generic
import os, sys, getopt
import matplotlib.pyplot as plt

# parse cmdline args in a C-like manner
from utils import get_closest_factors

try:
    opts, args = getopt.getopt(sys.argv[1:], "vd:i:o:n:t:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

verbose = False
dir_name = ""
trajs_dir_name = "ttrajs"
output_name = "sample_tic"
dim = 0

for o, a in opts:
    if o == "-v":
        verbose = True # do we even need this??
    elif o == "-d":
        dir_name = a
    elif o == "-o":
        output_name = a
    elif o == "-n":
        dim = int(a)
    elif o == "-t":
        trajs_dir_name = a

if dir_name != "":
    try:
        os.mkdir(dir_name)
    except:
        pass # ignore err due to existing dir
    os.chdir(dir_name)

meta, ttrajs = load_trajs(trajs_dir_name)
n_trajs = len(ttrajs.keys())
nx, ny = get_closest_factors(n_trajs)
fig1, axes = plt.subplots(nx, ny, figsize=(5 * nx, 3 * ny), sharex=True, sharey=True, squeeze=True)
for i, k in enumerate(ttrajs.keys()):
    traj = ttrajs[k]
    axes[i//ny][i % ny].plot(traj[:,dim])
    axes[i//ny][i % ny].set_ylim([-5,5]) # should extract this from data
fig1.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.xlabel("Time (ns)", fontsize=16)
plt.ylabel(f"tIC {dim}", fontsize=16)
plt.title(f"Projection of trajectories onto tIC {dim}", fontsize=16)
fig1.savefig(f"{output_name}_dim_{dim}.pdf")
