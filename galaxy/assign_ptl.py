#!/usr/bin/env python3
"""
Creates 2048^3 grids using the subhalo catalogue, separating them into red and blue populations.
"""

# import statements
import numpy as np
import h5py as hp
import sys
# import MAS_library as masl
import library_hicc.color as lhicc
import illustris_python as il
from library_hicc.printer import Printer

# reading command line inputs
SNAPSHOT = int(sys.argv[1])
BOX = int(sys.argv[2])


# defining needed paths
HOME = '/lustre/cosinga/tng%d/'%BOX
SAVE = '/lustre/cosinga/HI-color/results/fields/snap_%03d/'%SNAPSHOT
LOG = '/lustre/cosinga/HI-color/hicc/logs/galaxy/'

# getting simulation defined constants
f= hp.File(HOME+'snapdir_%03d/snap_%03d.0.hdf5'%(SNAPSHOT,SNAPSHOT), 'r')
head = dict(f['Header'].attrs)

LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
SCALE = head['Time'] # scale factor
BOXSIZE = head['BoxSize']/1e3 * SCALE #Mpc/h
REDSHIFT = head['Redshift']
DMPTL = head['MassTable'][1]*1e10/LITTLE_H

flds = ['SubhaloPos', 'SubhaloStellarPhotometrics', 'SubhaloMassType', 'SubhaloVel']
sub = il.groupcat.loadSubhalos(HOME,SNAPSHOT, fields=flds)
mass = sub[flds[2]][:]*1e10/LITTLE_H # solar masses
gr = sub[flds[1]][:,4] - sub[flds[1]][:,5]
pos = sub[flds[0]][:]/1e3 * SCALE # Mpc/h, 52 MB
vel = sub[flds[3]][:] # km/s, 52 MB
runs = ['low','mid','high']
w = hp.File(HOME + 'postprocessing/color_ptlID%d_%03d.hdf5'%(BOX,SNAPSHOT),'w')

for RUN in runs:
    resolved_mask = lhicc.is_resolved_stmass(mass[:,4])
    red_mask = lhicc.is_red_nelson(gr,mass[:,4],RUN)
    blue_mask = np.invert(red_mask)
    red_mask *= resolved_mask # removes unresolved true values
    blue_mask *= resolved_mask
    ptltypes = [0,1,4,5]
    blue_ptls = []
    red_ptls = []

    for i in range(len(blue_mask)):
        for p in ptltypes:
            ptl = il.snapshot.loadSubhalo(HOME, SNAPSHOT, i, p, fields='ParticleIDs')
            ptl = ptl[:]
            if blue_mask[i] == 1: # this subhalo is blue, label constituent ptls as blue
                blue_ptls.extend(list(ptl))
            elif red_mask[i] == 1:
                red_ptls.extend(list(ptl))
            
    w.create_dataset("blue_ptl_%s"%RUN,data=blue_ptls, compression="gzip", compression_opts=9)
    w.create_dataset("red_ptl_%s"%RUN, data=red_ptls, compression="gzip", compression_opts=9)
w.close()