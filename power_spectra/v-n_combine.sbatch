#!/bin/bash
# this file combines the v-n chunk files' into smaller chunks

#SBATCH --share
#SBATCH --job-name=v-n_combine
#SBATCH --time=1-5:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --requeue
#SBATCH --output=./outlogs/v-n_combine%a.dat
#SBATCH --error=./errors/v-n_combine%a.err

. ~/.profile
START=$SLURM_ARRAY_TASK_ID
END=$(($SLURM_ARRAY_TASK_ID+20))
python3 /lustre/cosinga/HI-color/hicc/vn/combine.py v-n $START $END $SNAP $BOX 0
python3 /lustre/cosinga/HI-color/hicc/vn/combine.py v-nrs $START $END $SNAP $BOX 0
