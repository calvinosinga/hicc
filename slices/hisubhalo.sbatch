#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hisubhalo_slice
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=15000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=hisubhalo_slice.dat

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/slices/slice.py hisubhalo 99 100
python3 /lustre/cosinga/HI-color/hicc/slices/slice.py hisubhalors 99 100
