#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hiptl-ptl_xpk
#SBATCH --time=2-10:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/hiptl-ptl_xpk.dat
#SBATCH --error=./errors/hiptl-ptl_xpk.err
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross $SNAP $BOX 1 $AXIS ptl hiptl
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross $SNAP $BOX 0 $AXIS ptlrs hiptlrs