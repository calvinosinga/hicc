#!/bin/bash

#SBATCH --share
#SBATCH --job-name=v-npk
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=./outlogs/v-npk.dat
#SBATCH --error=./errors/v-npk.err

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS v-n
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto $SNAP $BOX 0 $AXIS v-nrs