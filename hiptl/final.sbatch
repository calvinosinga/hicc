#!/bin/bash
# this file combines the hiptl chunk files' into the final field

#SBATCH --share
#SBATCH --job-name=hiptl_final
#SBATCH --output=logs/hiptl_final.log
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=ALL
#SBATCH --requeue

. ~/.profile
python3 /lustre/cosinga/hicc/vn/combine.py hiptl 0 448 99 100 1
