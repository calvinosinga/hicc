#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hisubhalo
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=hisubhalo.dat

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/hisubhalo/hisubhalo.py 99 100 1200 0
