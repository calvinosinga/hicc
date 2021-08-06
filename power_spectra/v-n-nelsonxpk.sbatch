#!/bin/bash

#SBATCH --share
#SBATCH --job-name=v-n-nelson_xpk
#SBATCH --time=1-10:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/v-n-nelson_xpk.dat
#SBATCH --error=./errors/v-n-nelson_xpk.err
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross $SNAP $BOX 1 $AXIS $RES nelson_mid v-n
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross $SNAP $BOX 0 $AXIS $RES nelsonrs_mid v-nrs
