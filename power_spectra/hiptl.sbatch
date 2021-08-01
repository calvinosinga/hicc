#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hiptl
#SBATCH --error=./errors/hiptl%a.err
#SBATCH --output=./outlogs/hiptl%a.log
#SBATCH --time=15:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --requeue

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/hiptl/hiptl.py $SLURM_ARRAY_TASK_ID $SNAP $BOX $RES $AXIS