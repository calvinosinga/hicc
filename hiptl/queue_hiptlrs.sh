#!/bin/bash
job1=$(sbatch hiptlrs.sbatch)
job2=$(sbatch --dependency=afterok:$job1 combine.sbatch)
sbatch --dependency=afterok:$job2 final.sbatch
