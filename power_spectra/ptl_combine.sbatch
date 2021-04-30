#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptl_combine
#SBATCH --output=logs/ptl_combine.log
#SBATCH --time=2-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=80000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/particle/combine.py ptl $SNAP $BOX $NUMFILES
python3 /lustre/cosinga/hicc/particle/combine.py ptlrs $SNAP $BOX $NUMFILES