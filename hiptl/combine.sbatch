#!/bin/bash
# this file combines the hiptl chunk files' into smaller chunks

#SBATCH --share
#SBATCH --array=0-440:20
#SBATCH --job-name=hiptl_combine
#SBATCH --output=logs/hiptl_combine%a.log
#SBATCH --time=21:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=72000
#SBATCH --mail-type=ALL
#SBATCH --requeue

. ~/.profile
START=$SLURM_ARRAY_TASK_ID
END=$(($SLURM_ARRAY_TASK_ID+20))
python3 /lustre/cosinga/HI-color/hicc/hiptl/combine.py hiptl $START $END 99 100 0
