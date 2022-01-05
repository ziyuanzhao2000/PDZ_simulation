#!/bin/bash
#SBATCH -c 1                # Number of cores (-c)
#SBATCH -t 0-01:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -serial_requeue
#SBATCH --mem=4G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o trajectory_to_snapshots.out
#SBATCH -e trajectory_to_snapshots.err
source ~/.bashrc
conda activate openmm
python trajectory_to_snapshots.py -d Dec15_liganded_PDZ -i production.h5 -o traj_snapshot -r reference.pdb -v