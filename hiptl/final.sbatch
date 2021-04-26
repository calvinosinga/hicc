#!/bin/bash
# this file combines the hiptl chunk files' into the final field

#SBATCH --share
#SBATCH --job-name=hiptl_final
#SBATCH --output=logs/hiptl_final.log
#SBATCH --time=2-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=80000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/hiptl/combine.py hiptl 0 448 99 100 1
#python3 /lustre/cosinga/hicc/hiptl/combine.py hiptl 0 448 67 100 1
