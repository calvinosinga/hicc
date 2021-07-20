#!/bin/bash

#SBATCH --share
#SBATCH --job-name=galaxy_slice
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=galaxy_slice.dat

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/slices/slice.py nelson_mid 99 100
python3 /lustre/cosinga/HI-color/hicc/slices/slice.py nelsonrs_mid 99 100
