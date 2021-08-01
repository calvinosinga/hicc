#!/bin/bash
#SBATCH --share
#SBATCH --job-name=hiptl_combine
#SBATCH --output=./outlogs/hiptl_combine%a.log
#SBATCH --error=./errors/hiptl_combine%a.err
#SBATCH --time=1-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --requeue

. ~/.profile
START=$SLURM_ARRAY_TASK_ID
END=$(($SLURM_ARRAY_TASK_ID+20))
python3 /lustre/cosinga/HI-color/hicc/hiptl/combine.py hiptl $START $END $SNAP $BOX 0
python3 /lustre/cosinga/HI-color/hicc/hiptl/combine.py hiptlrs $START $END $SNAP $BOX 0