#!/bin/bash
# this file combines the hiptl chunk files' into the final field

#SBATCH --share
#SBATCH --job-name=hiptlres_final
#SBATCH --output=hiptlres_final.log
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=100000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/hiptl/res_tests/combine.py hiptl 0 448 99 100 1
