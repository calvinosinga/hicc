#!/bin/bash

#SBATCH --share
#SBATCH --array=0-447
#SBATCH --job-name=hiptl
#SBATCH --error=./logs/hiptl%a.err
#SBATCH --output=./logs/hiptl%a.log
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=ALL
#SBATCH --requeue

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/hiptl/hiptl.py $SLURM_ARRAY_TASK_ID 99 100 1200 0
