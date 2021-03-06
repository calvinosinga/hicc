#!/usr/bin/env python3
"""
This file takes a few halos and subhalos of interest at redshift 0 and plots
the distribution of HI within them using the hiptl catalogue. The halos are
picked to have masses near 10^12 solar masses, and then are separated based
on whether they have a red or blue galaxy at its center.
"""

# import statements
import numpy as np
import h5py as hp
import illustris_python as il
import library_hicc.color as hcc
import matplotlib.pyplot as plt

# command-line inputs
CHUNK = sys.argv[1]
SNAPSHOT = sys.argv[2]
BOX = sys.argv[3]

# defining paths
HIPATH = '/lustre/diemer/illustris/hih2/'  # where the hiptl files are saved
TNG = '/lustre/cosinga/tng%s/'%BOX # where the ptl files are saved

# loading needed data
grp = il.loadHalos(TNG, int(SNAPSHOT), fields=['GroupFirstSub','Group_M_Mean200', 'GroupCM','Group_R_Crit200'])
sh = il.loadSubhalos(TNG, int(SNAPSHOT), fields=['SubhaloMassType','SubhaloCM','SubhaloStellarPhotometrics','SubhaloGrNr'])
head = il.loadHeader(TNG, int(SNAPSHOT))

# defining constants
LITTLE_H = head['HubbleParam']
BOXSIZE = head['BoxSize']/1e3 #Mpc/h

# getting the halos with masses near 10^12
target = 12
dist_to_target = np.abs(np.log10(grp['Group_M_Mean200'][:]*1e10/LITTLE_H) - target) # solar masses
targ_idx = dist_to_target < 0.5 # gets the halos that are within 11.5-12.5 range
print('We have found %d halos near the target...'%np.sum(targ_idx))

# getting the photometric data of the central and separating into blue/red
centrals_col = sh['SubhaloStellarPhotometrics'][grp['GroupFirstSub'][targ_idx], :]
centrals_stmass = sh['SubhaloMassType'][grp['GroupFirstSub'][targ_idx], 4] *1e10/LITTLE_H # solar masses
gr = centrals_col[:, 4] - centrals_col[:, 5]
red_idx = hcc.is_red_nelson(gr, centrals_stmass, 'mid')
blue_idx = ~red_idx
print('we have %d blues and %d reds found...'%(np.sum(blue_idx),np.sum(red_idx)))

# getting the group numbers of the halos with the first red and first blue centrals


# getting the ptl and hiptl data for those to compute the HI profiles

#
