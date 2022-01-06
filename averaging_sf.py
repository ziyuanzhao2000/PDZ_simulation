import reciprocalspaceship as rs
import numpy as np
import pandas as pd
import sys, os

# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "vDd:i:c:n")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

chain_id = 0
n_frames = 1
input_name = "traj_snapshot"
dirname = ""
verbose = False
dryrun = False

for o, a in opts:
    if o == "-v":
        verbose = True
    elif o == "-D":
        dryrun = True
    elif o == "-i":
        base_path = a
    elif o == "-d":
        dirname = a
    elif o == "-c":
        chain_id = int(a)
    elif o == "-n":
        n_frames = int(a)


if dirname != "":
    try:
        os.mkdir(dirname)
    except:
        pass # ignore err due to existing dir
    os.chdir(dirname)


try:
    dataset = rs.read_mtz(f"{input_name}_0_{chain_id}.mtz")
    n_reflections = dataset.shape[0]
except:
    print("Can't load the mtz file for the designated chain at frame 0!")
    sys.exit(1)

complex_reflections = np.zeros(n_reflections, dtype='complex128')

# sum up complex reflections and then take average
for frame in range(n_frames):
    dataset = rs.read_mtz(f"{input_name}_{frame}_{chain_id}.mtz")
    complex_reflections += np.array([amp*np.exp(np.pi*phase/180 * 1j) for [amp, phase] in dataset.to_numpy()])
    if verbose:
        print(f"Finished processing frame {frame}")
    if dryrun:
        sys.exit(1)
complex_reflections /= n_frames

# convert back to polar form and bind the composite indices
avg_dataset[:] = np.stack([np.abs(complex_reflections), np.angle(complex_reflections) / np.pi * 180]).T
avg_dataset.infer_mtz_dtypes(inplace = True)
avg_dataset.write_mtz(f"{input_name}_avg_{chain_id}.mtz")
if verbose:
    print("Wrote average structural factor to the .mtz file!")
