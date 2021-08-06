#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nelson_dust
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=nelson_dust.dat

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/galaxy/nelson_dust.py 99 100 'mid' 800 0
