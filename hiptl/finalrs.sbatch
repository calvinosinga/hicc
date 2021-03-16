#!/bin/bash
# this file combines the hiptl chunk files' into the final field

#SBATCH --share
#SBATCH --job-name=rshiptl_final
#SBATCH --output=logs/rshiptl_final.log
#SBATCH --time=2-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=80000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/hiptl/combine.py hiptlrs 0 448 99 100 1
