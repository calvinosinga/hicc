#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptl_final
#SBATCH --output=logs/ptl_final.log
#SBATCH --time=2-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/particle/combine.py ptl 0 448 99 100 1
