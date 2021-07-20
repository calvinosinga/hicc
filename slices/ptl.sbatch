#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptl_slice
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=ptl_slice.dat

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/slices/slice.py ptl 99 100
python3 /lustre/cosinga/HI-color/hicc/slices/slice.py ptlrs 99 100