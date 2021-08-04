#!/bin/bash
# this file combines the hiptl chunk files' into the final field

#SBATCH --share
#SBATCH --job-name=vnrs_final
#SBATCH --output=logs/vnrs_final.log
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=72000
#SBATCH --mail-type=ALL
#SBATCH --requeue

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/vn/combine.py v-nrs 0 448 99 100 1
