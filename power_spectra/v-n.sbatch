#!/bin/bash

#SBATCH --share
#SBATCH --array=0-447
#SBATCH --job-name=v-n
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --output=./outlogs/v-n%a.dat
#SBATCH --error=./errors/v-n%a.err

python3 /lustre/cosinga/HI-color/hicc/vn/vn.py $SLURM_ARRAY_TASK_ID $SNAP $BOX $RES $AXIS