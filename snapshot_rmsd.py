import numpy as np
import sys, os, getopt
import mdtraj as md
import matplotlib.pyplot as plt

# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "vDd:r:N:t:o:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

verbose = False
dryrun = False
dirname = ""
ref_name = ""
target_name = ""
output_name = "snapshots_rmsd"
n_frames = 100
delta_t = 1 # ns

for o, a in opts:
    if o == "-v":
        verbose = True
    elif o == "-D":
        dryrun = True
    elif o == "-d":
        dirname = a
    elif o == "-r":
        ref_name = a
    elif o == "-i":
        target_name = a
    elif o == "-N":
        n_frames = a
    elif o == "-t":
        delta_t = float(a)
    elif o == "-o":
        output_name = a

if dirname != "":
    try:
        os.mkdir(dirname)
    except:
        pass # ignore err due to existing dir
    os.chdir(dirname)

# load reference trajectory and topology from first frame
ref_traj = md.load(f"{ref_name}.pdb")
ref_top = ref_traj[0].topology
unit_cell_rmsd = []
for i in range(n_frames):
    target_traj = md.load(f"{target_name}_{i}.pdb")
    target_top = target_traj[0].topology
    unit_cell_rmsd.append(md.rmsd(target_traj, ref_traj, 0, ref_top.select("is_backbone")))

unit_cell_rmsd = np.array(unit_cell_rmsd)
plt.figure(figsize=(18, 12))
plt.rc('font', size=16)
plt.title("PDZ domain, unit cell RMSD")
plt.xlabel("Time (ns)")
plt.xticks(np.arange(0, n_frames, n_frames // 10), np.arange(0, n_frames, n_frames // 10) * delta_t)
plt.ylabel("RMSD (nm)")
plt.plot(unit_cell_rmsd, label = "unit cell 1")
plt.legend(loc = 'upper right')
plt.savefig(f'{output_name}.pdf')

