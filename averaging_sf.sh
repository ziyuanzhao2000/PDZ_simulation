#!/bin/bash
#SBATCH -c 1                # Number of cores (-c)
#SBATCH -t 0-01:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -serial_requeue
#SBATCH --mem=4G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o averaging_sf.out
#SBATCH -e averaging_sf.err
source ~/.bashrc
conda activate openmm # need to have reciprocalspaceship lib installed
python averaging_sf.py -d Dec15_liganded_PDZ -i traj_snapshot.h5 -c 0 -n 100 -v