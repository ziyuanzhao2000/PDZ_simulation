import mdtraj
import os
import sys, getopt

frame_skip = 10 # sample trajectory every 'frame_skip' frames
n_chains = 48 # 2 * 24

# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "vDcd:i:r:o:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

dirname = ""
verbose = False
dryrun = False
chainwise = False
inputname = "production.h5"
refname = "reference.pdb"
outputname = "liganded_pdz_traj_"


for o, a in opts:
    if o == "-v":
        verbose = True
    elif o == "-D":
        dryrun = True
    elif o == "-i":
        inputname = a
    elif o == "-o":
        outputname = a
    elif o == "-r":
        refname = a
    elif o == "-d":
        dirname = a
    elif o == "-c":
        chainwise = True

if dirname != "":
    try:
        os.mkdir(dirname)
    except:
        pass # ignore err due to existing dir
    os.chdir(dirname)

# Loads the .h5 trajectory file from the post-squeeze production run using OpenMM
target_traj = mdtraj.load(inputname)
top = target_traj.topology

# slice
target_traj = target_traj[::frame_skip]
ref_traj = mdtraj.load(refname)

if verbose:
    print("Loaded and sliced trajectory")
if dryrun:
    sys.exit(1)

for i, frame in enumerate(target_traj):
    if chainwise:
        for chain_id in range(0, n_chains, 2):
            chain = frame.atom_slice(top.select("chainid " + str(chain_id) + " or chainid " + str(chain_id + 1) + " and protein"))
            chain.superpose(ref_traj) # align using all atoms, water, ions, and protein included
            chain.save_pdb(f"{outputname}_{i}_{(chain_id / 2):d}.pdb")
    else:
        frame.superpose(ref_traj, atom_indices = top.select("is_backbone"), ref_atom_indices = top.select("is_backbone")).save_pdb(f"{outputname}_{i}.pdb")
    print(f"Finished outputting snapshots for frame {i * frame_skip}")