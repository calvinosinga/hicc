#!/bin/bash
# TODO: some way for rerunning individual jobs that fail in array in slurm
# TODO: error logs for array into one file - to find individual jobs that fail
# TODO: analysis step for paco, still need redshift space stuff from paco
# create directories to store output
mkdir outlogs
mkdir errors

rm ./errors/*.err
rm ./outlogs/*.dat
# check to make sure input is correct
if [ -z "$1" ]
then
    echo "$1 is not provided: needs snapshot 99,67,50,33"
    exit 125
else
    echo "snapshot: $1"
    
fi
if [ -z "$2" ]
then
    echo "$2 is not provided: needs box length: 100,300,50"
    exit 125
else
    echo "box length: $2"
fi
if [ -z "$3" ]
then
    echo "$3 is not provided: needs axis for line-of-sight - either 0,1,2 although typically 0"
    exit 125
else
    echo "axis: $3"
fi
if [ -z "$4" ]
then
    echo "command-line arg needed: number of files: 448-TNG100, 680-TNG50, 600-TNG300"
    exit 125
else
    echo "numfiles: $4"
fi
if [ -z "$5" ]
then
    echo "command-line arg needed: grid resolution"
    exit 125
else
    echo "resolution of grid: $5"
fi

ARRAYNO=$(($4-1))
COMBINE_ARRAY=$(($4-$4%20))

GRIDMEM=$(($5*$5*$5*4/1000000+10000))
PKMEM=$(($GRIDMEM*2))
XPKMEM=$(($PKMEM*2))

# submit hisubhalo jobs
hisubgrid=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --mem-per-cpu=$GRIDMEM hisubhalo.sbatch)
hisubgrid="${hisubgrid##* }"

sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$hisubgrid --mem-per-cpu=$PKMEM hisubhalopk.sbatch

# submit galaxy jobs
galgrid=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --mem-per-cpu=$GRIDMEM galaxy.sbatch)
galgrid="${galgrid##* }"

sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$galgrid --mem-per-cpu=$PKMEM galaxypk.sbatch

# nelson-nelson xpk
sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$galgrid --mem-per-cpu=$PKMEM nelson-nelsonxpk.sbatch

# hisubhalo-nelson xpk
sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$galgrid:$hisubgrid --mem-per-cpu=$XPKMEM hisubhalo-nelsonxpk.sbatch

# calculate pk for paco (need to change the memory)
# sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --mem-per-cpu=$PKMEM pacopk.sbatch

# paco-nelson xpk
# sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$galgrid --mem-per-cpu=$XPKMEM paco-nelsonxpk.sbatch

# submit hiptl jobs
hiptlgrid=$(sbatch --array=0-$ARRAYNO --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --mem-per-cpu=$GRIDMEM hiptl.sbatch)
hiptlgrid="${hiptlgrid##* }"

#SBATCH --array=0-460:20
hiptlcomb=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,NUMFILES=$4 --dependency=afterok:$hiptlgrid --array=0-$COMBINE_ARRAY:20 hiptl_combine.sbatch)
hiptlcomb="${hiptlcomb##* }"

hiptlfinal=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,NUMFILES=$4 --dependency=afterok:$hiptlcomb hiptl_final.sbatch)
hiptlfinal="${hiptlfinal##* }"

sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$hiptlfinal --mem-per-cpu=$PKMEM hiptlpk.sbatch

# submit ptl jobs
ptlgrid=$(sbatch --array=0-$ARRAYNO --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,RES=$5 --mem-per-cpu=$PKMEM ptl.sbatch)
ptlgrid="${ptlgrid##* }"

ptlcomb=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,NUMFILES=$4 --mem-per-cpu=$PKMEM --dependency=afterok:$ptlgrid --array=0-$COMBINE_ARRAY:20 ptl_combine.sbatch)
ptlcomb="${ptlcomb##* }"

ptlfinal=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,NUMFILES=$4 --mem-per-cpu=$PKMEM --dependency=afterok:$ptlcomb ptl_final.sbatch)
ptlfinal="${ptlfinal##* }"

sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$ptlfinal --mem-per-cpu=$PKMEM ptlpk.sbatch

# hiptl-nelson xpk
sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$galgrid:$hiptlcomb --mem-per-cpu=$XPKMEM hiptl-nelsonxpk.sbatch

# hiptl-ptl xpk
sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$hiptlcomb:$ptlcomb --mem-per-cpu=$XPKMEM hiptl-ptlxpk.sbatch


# paco-ptl xpk
# sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$ptlcomb --mem-per-cpu=$XPKMEM paco-ptlxpk.sbatch

# hisubhalo-ptl xpk
sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$hisubgrid:$ptlcomb --mem-per-cpu=$XPKMEM hisubhalo-ptlxpk.sbatch
# hydrotools visualization step?

# make gr-stmass plots
python /lustre/cosinga/hicc/color-stmass/gr-stmass.py $1 $2 &
