#!/bin/bash

#SBATCH --share
#SBATCH --job-name=vn_trial
#SBATCH --error=vn_trial.err
#SBATCH --output=vn_trial.log
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=27000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/vn/vn.py 12 99 100 1600 0