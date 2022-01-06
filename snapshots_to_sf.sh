#!/bin/bash
#SBATCH -c 1                # Number of cores (-c)
#SBATCH -t 0-01:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -serial_requeue
#SBATCH --mem=4G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o snapshots_to_sf.out
#SBATCH -e snapshots_to_sf.err
source ~/.bashrc
source /n/holylfs/LABS/hekstra_lab/garden/phenix/phenix-1.16-3549/phenix_env.sh
python snapshots_to_sf.py -d Dec15_liganded_PDZ -c 0 -f ${SLURM_ARRAY_TASK_ID} -v