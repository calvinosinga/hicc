#!/usr/bin/env python3
"""
Creates the 2048^3 grid for the HI/H2 catalogue assigned on a per
subhalo basis.
"""
# import statements
import numpy as np
import h5py as hp
import sys
print("about to import MAS_library...")
import MAS_library as masl
print("imported MAS_library")
from library_hicc.models import get_hisubhalo_models
import illustris_python as il

# reading command line inputs
SNAPSHOT = int(sys.argv[1])
BOX = int(sys.argv[2])

# defining needed paths
HOME = '/lustre/cosinga/tng%d/'%BOX
SAVE = '/lustre/cosinga/final_fields/'

# assigning author-defined constants (not expected to change)
GRID = (2048,2048,2048)
MAS = 'CIC'

# input data
f = hp.File(HOME+'/groups_%03d/hih2_galaxy_%03d.hdf5'%(SNAPSHOT,SNAPSHOT),'r')
pos = il.groupcat.loadSubhalos(HOME, SNAPSHOT, fields=['SubhaloCM'])
head = il.groupcat.loadHeader(HOME,SNAPSHOT)
pos = pos/1e3 # Mpc/h
models = get_hisubhalo_models()
print("the models used are: "+str(models)+'\n')

# output file
w = hp.File(SAVE+'hisubhalo%d_%03d.final.hdf5'%(BOX,SNAPSHOT), 'w')

# getting simulation-defined constants
BOXSIZE = head['BoxSize']/1e3 # Mpc/h

# loop over 9 models for HI
for m in models:
    print("computing the grid for %s"%m)
    field = np.zeros(GRID, dtype=np.float32)
    mass = f[m][:] # already in solar masses
    mass = mass.astype(np.float32)
    pos = pos.astype(np.float32)
    print("the sum is %f"%np.sum(mass))
    masl.MA(pos, field, BOXSIZE, MAS, mass)
    w.create_dataset(m, data=field, compression="gzip", compression_opts=9)

w.close()
