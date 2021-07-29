#!/bin/bash
# this file combines the hiptl chunk files' into smaller chunks

#SBATCH --share
#SBATCH --array=0-440:20
#SBATCH --job-name=v-nrs_combine
#SBATCH --output=logs/v-nrs_combine%a.log
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-type=ALL
#SBATCH --requeue

. ~/.profile
START=$SLURM_ARRAY_TASK_ID
END=$(($SLURM_ARRAY_TASK_ID+20))
python3 /lustre/cosinga/HI-color/hicc/vn/combine.py v-nrs $START $END 99 100 0