#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptl_combine%a
#SBATCH --output=./outlogs/ptl_combine%a.dat
#SBATCH --error=./errors/ptl_combine%a.err
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --requeue

. ~/.profile
START=$SLURM_ARRAY_TASK_ID
END=$(($SLURM_ARRAY_TASK_ID+20))
python3 /lustre/cosinga/HI-color/hicc/particle/combine.py ptl $START $END $SNAP $BOX 0
python3 /lustre/cosinga/HI-color/hicc/particle/combine.py ptlrs $START $END $SNAP $BOX 0
