#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hiptl-nelson_xpk
#SBATCH --time=2-10:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/hiptl-nelson_xpk.dat
#SBATCH --error=./errors/hiptl-nelson_xpk.err
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross $SNAP $BOX 1 $AXIS $RES nelson_mid hiptl
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross $SNAP $BOX 0 $AXIS $RES nelsonrs_mid hiptlrs