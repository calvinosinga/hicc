#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptl
#SBATCH --error=./errors/ptl%a.err
#SBATCH --output=./outlogs/ptl%a.log
#SBATCH --time=15:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/particle/particle.py $SLURM_ARRAY_TASK_ID $SNAP $BOX $RES -1
python3 /lustre/cosinga/hicc/particle/particle.py $SLURM_ARRAY_TASK_ID $SNAP $BOX $RES $AXIS
