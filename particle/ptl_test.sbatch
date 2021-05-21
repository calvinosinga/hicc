#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptlrs_test
#SBATCH --error=./logs/ptlrs_test.err
#SBATCH --output=./logs/ptlrs_test.log
#SBATCH --time=2-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=80000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/particle/particle.py 124 99 100 0
