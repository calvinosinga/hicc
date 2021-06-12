#!/bin/bash

#SBATCH --share
#SBATCH --job-name=pacopk
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=./outlogs/pacopk.dat
#SBATCH --error=./errors/pacopk.err

. ~/.profile
python3 /lustre/cosinga/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS paco