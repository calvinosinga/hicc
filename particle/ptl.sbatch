#!/bin/bash

#SBATCH --share
#SBATCH --array=0-447
#SBATCH --job-name=ptl
#SBATCH --error=./logs/ptl%a.err
#SBATCH --output=./logs/ptl%a.log
#SBATCH --time=2-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=5000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/particle/particle.py $SLURM_ARRAY_TASK_ID 99 100 800 0
