#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hiptl_slice
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=hiptl_slice.dat

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/slices/slice.py hiptl 99 100
python3 /lustre/cosinga/HI-color/hicc/slices/slice.py hiptlrs 99 100