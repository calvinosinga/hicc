#!/bin/bash

#SBATCH --share
#SBATCH --job-name=assign_ptl
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=assign_ptl.dat

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/galaxy/assign_ptl.py 99 100
