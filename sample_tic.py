# This script reads in the tICA-transformed trajectory stored in the ttrajs folder and then plot the time trace
# along each direction. But how to interpret this??

import numpy as np
from msmbuilder.io import load_trajs, load_generic
import os, sys, getopt
import matplotlib.pyplot as plt

# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "vd:i:o:n:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

verbose = False
dir_name = ""
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

if dir_name != "":
    try:
        os.mkdir(dir_name)
    except:
        pass # ignore err due to existing dir
    os.chdir(dir_name)

meta, ttrajs = load_trajs('ttrajs')
n_trajs = len(ttrajs.keys())

fig1 = plt.figure(figsize=(8 * 4, 6 * n_trajs))
for i, k in enumerate(ttrajs.keys()):
    traj = ttrajs[k]
    ax = fig1.add_subplot(n_trajs // 4, 4, i+1)
    ax.plot(traj[:,dim])
    ax.set_ylim([-3,3]) # should extract this from data
fig1.savefig(f"{output_name}_dim_{dim}.pdf")
