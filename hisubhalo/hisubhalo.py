#!/usr/bin/env python3
"""
Creates the 2048^3 grid for the HI/H2 catalogue assigned on a per
subhalo basis.
"""
# import statements
import numpy as np
import h5py as hp
import sys
from library_hicc.mas import CICW
from library_hicc.models import get_hisubhalo_models
import illustris_python as il
import redshift_space_library as rsl

# reading command line inputs, a third gives axis for redshift space
SNAPSHOT = int(sys.argv[1])
BOX = int(sys.argv[2])
RES = int(sys.argv[3]) # resolution of the grid
AXIS = int(sys.argv[4]) # is -1 if not in redshift space
IN_RS_SPACE = AXIS == -1


# defining needed paths
HOME = '/lustre/cosinga/tng%d/'%BOX
SAVE = '/lustre/cosinga/final_fields/'

# assigning author-defined constants (not expected to change)
GRID = (RES,RES,RES)

# getting simulation-defined constants
head = il.groupcat.loadHeader(HOME,SNAPSHOT)
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
BOXSIZE = head['BoxSize']/1e3 #Mpc/h
REDSHIFT = head['Redshift']
SCALE = head['Time'] # scale factor

print('the boxsize is %f'%BOXSIZE)

# input data
f = hp.File(HOME+'/groups_%03d/hih2_galaxy_%03d.hdf5'%(SNAPSHOT,SNAPSHOT),'r')
ids = f['id_subhalo'][:] # used to idx into the subhalo catalog
ids = ids.astype(np.int32)
sub = il.groupcat.loadSubhalos(HOME, SNAPSHOT, fields=['SubhaloPos','SubhaloVel'])
pos = sub['SubhaloPos'][ids]/1e3 * SCALE # Mpc/h
vel = sub['SubhaloVel'][ids] # km/s
models = get_hisubhalo_models()
print("the models used are: "+str(models)+'\n')
del head, sub

# checking to make sure that the positions match up with previous version


# if we are in redshift space, shift positions using velocity
# then make output file, changing the name depending on if its in real space or redshift space
if IN_RS_SPACE:
    rsl.pos_redshift_space(pos, vel, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)
    w = hp.File(SAVE+"hisubhalors%d_%03d.final.hdf5"%(BOX,SNAPSHOT), 'w')
else:
    w = hp.File(SAVE+'hisubhalo%d_%03d.final.hdf5'%(BOX,SNAPSHOT), 'w')

# loop over 9 models for HI
for m in models:
    print("computing the grid for %s"%m)
    field = np.zeros(GRID, dtype=np.float32) # ~32GB
    mass = f[m][:] # already in solar masses
    mass = mass.astype(np.float32)
    pos = pos.astype(np.float32)
    print("the average is %f solar masses"%np.mean(mass))
    CICW(pos, field, BOXSIZE, mass)
    w.create_dataset(m, data=field, compression="gzip", compression_opts=9)

w.close()
f.close()
