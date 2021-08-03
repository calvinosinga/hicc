#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hisubhalo_pk
#SBATCH --error=./errors/hisubhalo_pk.err
#SBATCH --time=2-18:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/hisubhalo_pk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS hisubhalo
printf "finished hisubhalo pk calculation \n"

python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto $SNAP $BOX 0 $AXIS hisubhalors
printf "finished hisubhalors pk calculation \n"

printf "now starting correlation calculations\n"

python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py auto $SNAP $BOX $AXIS hisubhalo
printf "finished hisubhalo correlation\n"

python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py auto $SNAP $BOX $AXIS hisubhalors
printf "finished hisubhalors correlation\n"

printf "now making slice plots\n"

python3 /lustre/cosinga/HI-color/hicc/slices/slice.py hisubhalo $SNAP $BOX
printf "made hisubhalo slice plot\n"

python3 /lustre/cosinga/HI-color/hicc/slices/slice.py hisubhalors $SNAP $BOX
printf "made hisubhalors slice plot\n"

printf "finished hisubhalo analysis!"
