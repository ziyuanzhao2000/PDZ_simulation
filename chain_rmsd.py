import mdtraj
import matplotlib.pyplot as pyplot
import os
import sys, getopt

try:
    opts, args = getopt.getopt(sys.argv[1:], "Dvd:i:o:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)


dirname = ""
verbose = False
dryrun = False
inputname = "iter"
outputname = "squeeze_analysis_plot.png"

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

target_traj = mdtraj.load('xxx.h5')
ref_traj = mdtraj.load('xxx')


# Calculate rmsd
top = ref_traj.topology

if dryrun:
    sys.exit(1)
    
for each chain in top.chains:
    chain_rmsd = mdtraj.rmsd(traj, reference, 0, chain)
if verbose:
    print("Calculated RMSD for all chains")
