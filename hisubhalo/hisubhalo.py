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
from library_hicc.redshift_space import pos_redshift_space
from library_hicc.printer import Printer

# reading command line inputs, a third gives axis for redshift space
SNAPSHOT = int(sys.argv[1])
BOX = int(sys.argv[2])
RES = int(sys.argv[3]) # resolution of the grid
AXIS = int(sys.argv[4]) # is -1 if not in redshift space


# defining needed paths
HOME = '/lustre/cosinga/tng%d/'%BOX
SAVE = '/lustre/cosinga/HI-color/results/fields/snap_%03d/'%SNAPSHOT
LOG = '/lustre/cosinga/HI-color/hicc/logs/hisubhalo/'

# assigning author-defined constants (not expected to change)
GRID = (RES,RES,RES)

# getting simulation-defined constants
head = il.groupcat.loadHeader(HOME,SNAPSHOT)
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
SCALE = head['Time'] # scale factor
BOXSIZE = head['BoxSize']/1e3 * SCALE #Mpc/h
REDSHIFT = head['Redshift']

# if we are in redshift space, shift positions using velocity
# then make output file, changing the name depending on if its in real space or redshift space
wrs = hp.File(SAVE+"hisubhalors%d_%03d.final.hdf5"%(BOX,SNAPSHOT), 'w')
w = hp.File(SAVE+'hisubhalo%d_%03d.final.hdf5'%(BOX,SNAPSHOT), 'w')
pnt = Printer(LOG+'hisubhalo%d_%03d.log'%(BOX,SNAPSHOT))
pnt.write('the boxsize is %f'%BOXSIZE)

# input data
f = hp.File(HOME+'/groups_%03d/hih2_galaxy_%03d.hdf5'%(SNAPSHOT,SNAPSHOT),'r')
ids = f['id_subhalo'][:] # used to idx into the subhalo catalog
ids = ids.astype(np.int32)
sub = il.groupcat.loadSubhalos(HOME, SNAPSHOT, fields=['SubhaloPos','SubhaloVel'])
pos = sub['SubhaloPos'][ids]/1e3 * SCALE # Mpc/h
vel = sub['SubhaloVel'][ids] # km/s
models = get_hisubhalo_models()
pnt.write("the models used are: "+str(models)+'\n')
rspos = pos_redshift_space(pos, vel, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)
del head, sub

# loop over 9 models for HI
for m in models:
    pnt.write("computing the grid for %s"%m)
    field = np.zeros(GRID, dtype=np.float32) # ~32GB
    mass = f[m][:] # already in solar masses
    mass = mass.astype(np.float32)
    pos = pos.astype(np.float32)

    pnt.writeTab("now assigning to the real-space grid")

    CICW(pos, field, BOXSIZE, mass)

    pnt.writeTab("now saving the grid to the file...")
    w.create_dataset(m, data=field, compression="gzip", compression_opts=9)

    pnt.writeTab("now creating new grid for redshift space...")
    field = np.zeros(GRID, dtype=np.float32) # ~32GB
    pnt.writeTab("assigning to the redshift-grid...")
    CICW(rspos, field, BOXSIZE, mass)
    pnt.writeTab("now saving to redshift file...")
    wrs.create_dataset(m,data=field, compression="gzip", compression_opts=9)


wrs.close()
w.close()
f.close()
