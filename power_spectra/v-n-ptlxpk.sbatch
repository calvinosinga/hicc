#!/bin/bash

#SBATCH --share
#SBATCH --job-name=v-n-ptl_xpk
#SBATCH --time=1-10:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/v-n-ptl_xpk.dat
#SBATCH --error=./errors/v-n-ptl_xpk.err
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross $SNAP $BOX 1 $AXIS ptl v-n
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross $SNAP $BOX 1 $AXIS ptlrs v-nrs