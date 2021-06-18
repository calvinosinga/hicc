#!/bin/bash

#SBATCH --share
#SBATCH --array=0-440:20
#SBATCH --job-name=hiptl_orgtest
#SBATCH --output=hiptl_orgtest%a.log
#SBATCH --time=1-10:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=50000
#SBATCH --mail-type=ALL

. ~/.profile
START=$SLURM_ARRAY_TASK_ID
END=$(($SLURM_ARRAY_TASK_ID+20))
python3 /lustre/cosinga/hicc/hiptl/org_test/hiptl.py $START $END 99 100 1600 0