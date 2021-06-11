#!/bin/bash
#SBATCH --share
#SBATCH --array=0-440:20
#SBATCH --job-name=hiptlrs_combine
#SBATCH --output=logs/hiptlrs_combine%a.log
#SBATCH --time=2-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=50000
#SBATCH --mail-type=ALL

. ~/.profile
START=$SLURM_ARRAY_TASK_ID
END=$(($SLURM_ARRAY_TASK_ID+20))
python3 /lustre/cosinga/hicc/hiptl/combine.py hiptlrs $START $END 99 100 0
