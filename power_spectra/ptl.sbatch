#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptl
#SBATCH --error=./logs/ptl%a.err
#SBATCH --output=./logs/ptl%a.log
#SBATCH --time=15:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=80000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/particle/ptl.py $SLURM_ARRAY_TASK_ID $SNAP $BOX
python3 /lustre/cosinga/hicc/particle/ptl.py $SLURM_ARRAY_TASK_ID $SNAP $BOX $AXIS