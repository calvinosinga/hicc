#!/bin/bash
# TODO: some way for rerunning individual jobs that fail in array in slurm
# TODO: error logs for array into one file - to find individual jobs that fail
# TODO: analysis step for paco, still need redshift space stuff from paco
# create directories to store output
mkdir outlogs
mkdir errors

rm ./errors/*.err
rm ./outlogs/*.dat

rm /lustre/cosinga/HI-color/hicc/logs/corr/*
rm /lustre/cosinga/HI-color/hicc/logs/galaxy/*
rm /lustre/cosinga/HI-color/hicc/logs/hiptl/*
rm /lustre/cosinga/HI-color/hicc/logs/hisubhalo/*
rm /lustre/cosinga/HI-color/hicc/logs/pk/*
rm /lustre/cosinga/HI-color/hicc/logs/ptl/*
rm /lustre/cosinga/HI-color/hicc/logs/slices/*
rm /lustre/cosinga/HI-color/hicc/logs/v-n/*

# check to make sure input is correct

printf "usage: bash calculate_pk.sh [snapshot] [boxsize] [axis] [numfiles] [grid_resolution]\n"
if [ -z "$1" ]
then
    printf "$1 is not provided: needs snapshot 99,67,50,33\n"
    exit 125
else
    printf "snapshot: $1\n"
    
fi
if [ -z "$2" ]
then
    printf "$2 is not provided: needs box length: 100,300,50\n"
    exit 125
else
    printf "box length: $2\n"
fi
if [ -z "$3" ]
then
    printf "$3 is not provided: needs axis for line-of-sight - either 0,1,2 although typically 0\n"
    exit 125
else
    printf "axis: $3\n"
fi
if [ -z "$4" ]
then
    printf "command-line arg needed: number of files: 448-TNG100, 680-TNG50, 600-TNG300 \n"
    exit 125
else
    printf "numfiles: $4\n"
fi
if [ -z "$5" ]
then
    printf "command-line arg needed: grid resolution\n"
    exit 125
else
    printf "resolution of grid: $5\n\n"
fi

ARRAYNO=$(($4-1))
COMBINE_ARRAY=$(($4-$4%20))

GRIDMEM=$(($5*$5*$5*4/1000000+10000))
COMBINEMEM=$(($GRIDMEM*2))
PKMEM=$(($GRIDMEM*2+15000))
XPKMEM=$(($PKMEM*2))

printf "The number of array jobs (should be one less than numfiles): $ARRAYNO\n"
printf "The number of combine arrays: $COMBINE_ARRAY\n\n"

printf "The memory for jobs that create a grid: $GRIDMEM\n"
printf "The memory for jobs that combine grids: $COMBINEMEM\n"
printf "The memory for jobs that calculate an auto power spectrum:$PKMEM\n"
printf "The memory for jobs that calculate a cross-power: $XPKMEM\n"

# submit hisubhalo jobs
hisubgrid=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --mem-per-cpu=$GRIDMEM hisubhalo.sbatch)
hisubgrid="${hisubgrid##* }"

sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --dependency=afterok:$hisubgrid --mem-per-cpu=$PKMEM hisubhalopk.sbatch

# submit galaxy jobs
galgrid=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --mem-per-cpu=$GRIDMEM galaxy.sbatch)
galgrid="${galgrid##* }"

sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --dependency=afterok:$galgrid --mem-per-cpu=$PKMEM galaxypk.sbatch

# nelson-nelson xpk
sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --dependency=afterok:$galgrid --mem-per-cpu=$PKMEM nelson-nelsonxpk.sbatch

# hisubhalo-nelson xpk
sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --dependency=afterok:$galgrid:$hisubgrid --mem-per-cpu=$XPKMEM hisubhalo-nelsonxpk.sbatch

# submit hiptl jobs
hiptlgrid=$(sbatch --array=0-$ARRAYNO --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --mem-per-cpu=$GRIDMEM hiptl.sbatch)
hiptlgrid="${hiptlgrid##* }"

hiptlcomb=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,NUMFILES=$4 --dependency=afterok:$hiptlgrid --array=0-$COMBINE_ARRAY:20 --mem-per-cpu=$COMBINEMEM hiptl_combine.sbatch)
hiptlcomb="${hiptlcomb##* }"

hiptlfinal=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,NUMFILES=$4 --dependency=afterok:$hiptlcomb --mem-per-cpu=$COMBINEMEM hiptl_final.sbatch)
hiptlfinal="${hiptlfinal##* }"

sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --dependency=afterok:$hiptlfinal --mem-per-cpu=$PKMEM hiptlpk.sbatch

# submit ptl jobs
ptlgrid=$(sbatch --array=0-$ARRAYNO --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --mem-per-cpu=$GRIDMEM ptl.sbatch)
ptlgrid="${ptlgrid##* }"

ptlcomb=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,NUMFILES=$4 --mem-per-cpu=$COMBINEMEM --dependency=afterok:$ptlgrid --array=0-$COMBINE_ARRAY:20 ptl_combine.sbatch)
ptlcomb="${ptlcomb##* }"

ptlfinal=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,NUMFILES=$4 --mem-per-cpu=$COMBINEMEM --dependency=afterok:$ptlcomb ptl_final.sbatch)
ptlfinal="${ptlfinal##* }"

sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --dependency=afterok:$ptlfinal --mem-per-cpu=$PKMEM ptlpk.sbatch

# submit v-n jobs
vngrid=$(sbatch --array=0-$ARRAYNO --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --mem-per-cpu=$GRIDMEM v-n.sbatch)
vngrid="${vngrid##* }"

vncomb=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,NUMFILES=$4 --mem-per-cpu=$COMBINEMEM --dependency=afterok:$vngrid --array=0-$COMBINE_ARRAY:20 v-n_combine.sbatch)
vncomb="${vncomb##* }"

vnfinal=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,NUMFILES=$4 --mem-per-cpu=$COMBINEMEM --dependency=afterok:$vncomb v-n_final.sbatch)
vnfinal="${vnfinal##* }"

sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --dependency=afterok:$vnfinal --mem-per-cpu=$PKMEM v-npk.sbatch

# hiptl-nelson xpk
sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --dependency=afterok:$galgrid:$hiptlfinal --mem-per-cpu=$XPKMEM hiptl-nelsonxpk.sbatch

# hiptl-ptl xpk
sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --dependency=afterok:$hiptlfinal:$ptlfinal --mem-per-cpu=$XPKMEM hiptl-ptlxpk.sbatch

# hisubhalo-ptl xpk
sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --dependency=afterok:$hisubgrid:$ptlfinal --mem-per-cpu=$XPKMEM hisubhalo-ptlxpk.sbatch

# v-n-ptl xpk
sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --dependency=afterok:$vnfinal:$ptlfinal --mem-per-cpu=$XPKMEM v-n-ptlxpk.sbatch

# v-n-galaxy xpk
sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --dependency=afterok:$vnfinal:$galgrid --mem-per-cpu=$XPKMEM v-n-nelsonxpk.sbatch

# make gr-stmass plots
python /lustre/cosinga/HI-color/hicc/color-stmass/gr-stmass.py $1 $2 $3 &
