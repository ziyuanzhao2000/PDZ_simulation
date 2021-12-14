#!/bin/bash
#SBATCH -c 1                # Number of cores (-c)
#SBATCH -t 0-00:10          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p test
#SBATCH --mem=100M           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o pdz_crystal_MD.out  
#SBATCH -e pdz_crystal_MD.err 

source ~/.bashrc
conda activate openmm
python squeeze.py -D -d Dec14
