#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH -t 7-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p gpu
#SBATCH --mem=4G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o pdz_solvated_sim_gpu.out
#SBATCH -e pdz_solvated_sim_gpu.err

source ~/.bashrc
conda activate openmm
python simple_run.py -d $1 -n $2 -t $3 -v
