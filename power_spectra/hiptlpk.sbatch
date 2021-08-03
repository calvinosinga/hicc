#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hiptl_pk
#SBATCH --error=./errors/hiptl_pk.err
#SBATCH --time=2-18:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/hiptl_pk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS hiptl
printf "finished hiptl pk calculation\n"

python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto $SNAP $BOX 0 $AXIS hiptlrs
printf "finished hiptlrs pk calculation\n"

printf "now calculating correlations\n"

python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py auto $SNAP $BOX $AXIS hiptl
printf "finished hiptl correlation\n"

python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py auto $SNAP $BOX $AXIS hiptlrs
printf "finished hiptlrs correlation\n"

printf "now starting to make slices...\n"

python3 /lustre/cosinga/HI-color/hicc/slices/slice.py hiptl $SNAP $BOX
printf "finished hiptl slice plot\n"

python3 /lustre/cosinga/HI-color/hicc/slices/slice.py hiptlrs $SNAP $BOX
printf "finished hiptlrs slice plot\n"

printf "finished hiptl output!"