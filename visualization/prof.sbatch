#!/bin/bash

#SBATCH --share
#SBATCH --job-name=prof
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30000
#SBATCH --account=astronomy
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=prof.dat

. ~/.profile
python3 /lustre/cosinga/hicc/visualization/ht_prof.py
