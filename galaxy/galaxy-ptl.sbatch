#!/bin/bash

#SBATCH --share
#SBATCH --job-name=galaxy-ptl
#SBATCH --time=1-10:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=galaxy-ptl.dat

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/galaxy/galaxy-ptl.py 99 100 'mid' 800 0
