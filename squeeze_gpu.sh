#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH -t 7-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p gpu
#SBATCH --mem=128G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o pdz_crystal_MD.out  
#SBATCH -e pdz_crystal_MD.err 

source ~/.bashrc
conda activate openmm
python squeeze.py -d Dec14
