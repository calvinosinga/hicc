#!/bin/bash

#SBATCH --share
#SBATCH --error=./errors/hisubhalo.err
#SBATCH --job-name=hisubhalo
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=./outlogs/hisubhalo.dat

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/hisubhalo/hisubhalo.py $SNAP $BOX $RES $AXIS
