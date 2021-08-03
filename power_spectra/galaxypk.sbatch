#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nelson_pk
#SBATCH --error=./errors/nelson_pk.err
#SBATCH --time=2-20:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/nelson_pk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS nelson_low
printf "finished auto pk calculation for nelson_low\n"

python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS nelson_high
printf "finished auto pk calculation for nelson_high\n"

python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS nelson_mid
printf "finished auto pk calculation for nelson_mid\n"

python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto $SNAP $BOX 0 $AXIS nelsonrs_mid
printf "finished auto pk calculation for nelson_mid in redshift-space\n"

printf "now starting correlation calculations\n"

python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py auto $SNAP $BOX $AXIS nelson_mid
printf "finished auto correlation calculation for nelson_mid\n"

python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py auto $SNAP $BOX $AXIS nelsonrs_mid
printf "finished auto correlation calculation for nelson_mid in redshift-space\n"

printf "now making slice plots...\n"

python3 /lustre/cosinga/HI-color/hicc/slices/slice.py auto nelson_mid $SNAP $BOX 
printf "made slice plots for nelson_mid\n"

python3 /lustre/cosinga/HI-color/hicc/slices/slice.py auto nelsonrs_mid $SNAP $BOX
printf "made slice plots for nelsonrs_mid\n"

printf "finished operations for galaxies!"