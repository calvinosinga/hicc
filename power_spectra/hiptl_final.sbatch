#!/bin/bash
# this file combines the hiptl chunk files' into the final field

#SBATCH --share
#SBATCH --job-name=hiptl_final
#SBATCH --output=./outlogs/hiptl_final.log
#SBATCH --error=./errors/hiptl_final.err
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/hiptl/combine.py hiptl 0 $NUMFILES $SNAP $BOX 1
python3 /lustre/cosinga/hicc/hiptl/combine.py hiptlrs 0 $NUMFILES $SNAP $BOX 1