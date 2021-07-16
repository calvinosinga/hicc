#!/bin/bash

#SBATCH --share
#SBATCH --array=0-447
#SBATCH --job-name=v-n
#SBATCH --output=./logs/v-n%a.log
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=5000
#SBATCH --mail-type=ALL
#SBATCH --requeue

python3 /lustre/cosinga/hicc/vn/vn.py $SLURM_ARRAY_TASK_ID 99 100 800 0
