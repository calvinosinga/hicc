#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hisubhalo-ptl_xpk
#SBATCH --time=2-10:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/hisubhalo-ptl_xpk.dat
#SBATCH --error=./errors/hisubhalo-ptl_xpk.err
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross $SNAP $BOX 1 $AXIS $RES hisubhalo ptl
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross $SNAP $BOX 0 $AXIS $RES hisubhalors ptlrs
