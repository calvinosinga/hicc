#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptl_combine
#SBATCH --output=./outlogs/ptl_combine.dat
#SBATCH --error=./errors/ptl_combine.err
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL

. ~/.profile
START=$SLURM_ARRAY_TASK_ID
END=$(($SLURM_ARRAY_TASK_ID+20))
python3 /lustre/cosinga/hicc/particle/combine.py ptl $START $END $SNAP $BOX 0
python3 /lustre/cosinga/hicc/particle/combine.py ptlrs $START $END $SNAP $BOX 1
