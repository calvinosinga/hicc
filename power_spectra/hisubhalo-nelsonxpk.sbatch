#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hisubhalo-nelson_xpk
#SBATCH --time=2-10:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/hisubhalo-nelson_xpk.dat
#SBATCH --error=./errors/hisubhalo-nelson_xpk.err
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/pk/pk.py cross $SNAP $BOX 1 $AXIS hisubhalo nelson_mid 0
python3 /lustre/cosinga/hicc/pk/pk.py cross $SNAP $BOX 0 $AXIS hisubhalo nelsonrs_mid 0