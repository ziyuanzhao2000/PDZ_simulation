#!/bin/bash
#SBATCH -c 1                # Number of cores (-c)
#SBATCH -t 0-01:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -serial_requeue
#SBATCH --mem=4G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o snapshots_to_sf.out
#SBATCH -e snapshots_to_sf.err
source ~/.bashrc
conda activate openmm
python snapshots_to_sf.py -d Dec15_liganded_PDZ -I ${SLURM_ARRAY_TASK_ID} -v