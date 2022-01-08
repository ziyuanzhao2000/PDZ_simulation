# This script works with slurm job batcher to
# convert input snapshots of MD trajectory to
# structural factors using phenix.fmodel

import sys, os, getopt

try:
    opts, args = getopt.getopt(sys.argv[1:], "Cvi:I:f:c:N:n:d:r:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

verbose = False
chainwise = False
input_name = "" # must supply input name!
dirname = ""
frame_id = 0
chain_id = 0
array_index = 0
n_frame = 100
n_chain = 24
resolution = 1.5

for o, a in opts:
    if o == "-v":
        verbose = True
    elif o == "-i":
        input_name = a
    elif o == "-f":
        frame_id = int(a)
    elif o == "-c":
        chain_id = int(a)
    elif o == "-I":
        array_index = int(a)
    elif o == "-N":
        n_frame = int(a)
    elif o == "-n":
        n_chain = int(a)
    elif o == "-d":
        dirname = a
    elif o == "-C":
        chainwise = True
    elif o == "-r":
        resolution = float(a)

if dirname != "":
    try:
        os.mkdir(dirname)
    except:
        pass # ignore err due to existing dir
    os.chdir(dirname)


if array_index > 0:
    frame_id = array_index // n_chain
    chain_id = array_index % n_chain

if frame_id < 0 or frame_id >= n_frame or chain_id < 0 or chain_id >= n_chain:
    print("Please input a correct frame and chain id!")
    sys.exit(1)

if verbose:
    print("Starting to call phenix.fmodel.")
if chainwise:
    os.system(f'phenix.fmodel {input_name}_{frame_id}_{chain_id}.pdb high_resolution={resolution}')
else:
    os.system(f'phenix.fmodel {input_name}_{frame_id}.pdb high_resolution={resolution}')

if verbose:
    print("Phenix.fmodel finished.")