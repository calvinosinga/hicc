#!/usr/bin/env python3
"""
Creates the 2048^3 grid for the HI/H2 catalogue assigned on a per
subhalo basis.
"""
# import statements
print("starting import statements")
import numpy as np
import h5py as hp
import sys
print("about to import MAS_library...")
import MAS_library as masl
print("imported MAS_library")
from library_hicc.models import get_hisubhalo_models
import illustris_python as il
import redshift_space_library as rsl

# reading command line inputs, a third gives axis for redshift space
SNAPSHOT = int(sys.argv[1])
BOX = int(sys.argv[2])
if len(sys.argv) > 3:
    AXIS = int(sys.argv[3])
    IN_RS_SPACE = True
else:
    IN_RS_SPACE = False

# defining needed paths
HOME = '/lustre/cosinga/tng%d/'%BOX
SAVE = '/lustre/cosinga/final_fields/'

# assigning author-defined constants (not expected to change)
GRID = (2048,2048,2048)
MAS = 'CIC'

# getting simulation-defined constants
head = il.groupcat.loadHeader(HOME,SNAPSHOT)
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
BOXSIZE = head['BoxSize']/1e3 #Mpc/h
REDSHIFT = head['Redshift']

# input data
f = hp.File(HOME+'/groups_%03d/hih2_galaxy_%03d.hdf5'%(SNAPSHOT,SNAPSHOT),'r')
sub = il.groupcat.loadSubhalos(HOME, SNAPSHOT, fields=['SubhaloCM','SubhaloVel'])
pos = sub['SubhaloCM']/1e3 # Mpc/h
vel = sub['SubhaloVel'] # km/s
models = get_hisubhalo_models()
print("the models used are: "+str(models)+'\n')
del head, sub

# if we are in redshift space, shift positions using velocity
if IN_RS_SPACE:
    rsl.pos_redshift_space(pos, vel, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)

# output file
w = hp.File(SAVE+'hisubhalo%d_%03d.final.hdf5'%(BOX,SNAPSHOT), 'w')

# loop over 9 models for HI
for m in models:
    print("computing the grid for %s"%m)
    field = np.zeros(GRID, dtype=np.float32) # ~32GB
    mass = f[m][:] # already in solar masses
    mass = mass.astype(np.float32)
    pos = pos.astype(np.float32)
    print("the sum is %f"%np.sum(mass))
    masl.MA(pos, field, BOXSIZE, MAS, mass)
    w.create_dataset(m, data=field, compression="gzip", compression_opts=9)

w.close()
f.close()