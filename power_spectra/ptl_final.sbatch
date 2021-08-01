#!/bin/bash
# this file combines the hiptl chunk files' into the final field

#SBATCH --share
#SBATCH --job-name=ptl_final
#SBATCH --output=./outlogs/ptl_final.log
#SBATCH --error=./errors/ptl_final.err
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/particle/combine.py ptl 0 $NUMFILES $SNAP $BOX 1
python3 /lustre/cosinga/HI-color/hicc/particle/combine.py ptlrs 0 $NUMFILES $SNAP $BOX 1