#!/bin/bash
# this file combines the hiptl chunk files' into smaller chunks

#SBATCH --share
#SBATCH --array=0-440:20
#SBATCH --job-name=v-n_combine
#SBATCH --output=logs/v-n_combine%a.log
#SBATCH --time=2-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=ALL
#SBATCH --requeue

. ~/.profile
START=$SLURM_ARRAY_TASK_ID
END=$(($SLURM_ARRAY_TASK_ID+20))
python3 /lustre/cosinga/hicc/vn/combine.py v-n $START $END 99 100 0