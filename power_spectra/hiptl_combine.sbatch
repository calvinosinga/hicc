#!/bin/bash
#SBATCH --share
#SBATCH --job-name=hiptl_combine
#SBATCH --output=./outlogs/hiptl_combine.log
#SBATCH --error=./errors/hiptl_combine.err
#SBATCH --time=2-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=80000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/hiptl/combine.py hiptl $SNAP $BOX $NUMFILES
python3 /lustre/cosinga/hicc/hiptl/combine.py hiptlrs $SNAP $BOX $NUMFILES