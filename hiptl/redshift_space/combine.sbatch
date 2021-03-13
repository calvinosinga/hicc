#!/bin/bash
# this file combines the hiptl chunk files' fields into smaller chunks
# to be combined again later

#SBATCH --share
#SBATCH --array=0-460:20
#SBATCH --job-name=hiptlrs_combine
#SBATCH --output=./logs/rshiptl_combine%a.log
#SBATCH --time=2-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=80000
#SBATCH --mail-type=ALL

. ~/.profile
START=$SLURM_ARRAY_TASK_ID
END=$(($SLURM_ARRAY_TASK_ID+20))
python3 /lustre/cosinga/hicc/hiptl/redshift_space/combine.py $START $END 99 100 0
