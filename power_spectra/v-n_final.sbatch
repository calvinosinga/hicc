#!/bin/bash
# this file combines the v-n chunk files' into smaller chunks

#SBATCH --share
#SBATCH --job-name=v-n_final
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --output=./outlogs/v-n_final.dat
#SBATCH --error=./errors/v-n_final.err

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/vn/combine.py v-n 0 $NUM_FILES $SNAP $BOX 1
python3 /lustre/cosinga/HI-color/hicc/vn/combine.py v-nrs 0 $NUM_FILES $SNAP $BOX 1
