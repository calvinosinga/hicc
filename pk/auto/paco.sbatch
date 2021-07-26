#!/bin/bash

#SBATCH --share
#SBATCH --job-name=pacopk
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=170000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=pacopk.dat

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto 99 100 1 0 paco
