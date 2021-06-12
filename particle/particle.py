#!/usr/bin/env python3
"""
This file takes each of the 448 subfiles for the particle subfiles and creates a 2048^3
grid for them.
"""
# import statements
import numpy as np
import h5py as hp
import sys
from library_hicc.mas import CICW
#import redshift_space_library as rsl
from library_hicc.redshift_space import pos_redshift_space

# reading command line inputs
CHUNK = int(sys.argv[1])
SNAPSHOT = int(sys.argv[2])
BOX = int(sys.argv[3])
RES = int(sys.argv[4]) # resolution of the grid
AXIS = int(sys.argv[5]) # if -1, not in redshift space
IN_RS_SPACE = not (AXIS == -1)

# defining needed paths
PTLPATH = '/lustre/cosinga/tng%d/snapdir_%03d/'%(BOX, SNAPSHOT) # where the ptl files are saved
OUTPATH = '/lustre/cosinga/hiptl_output/' # where to save the output
HIPATH = '/lustre/diemer/illustris/hih2/'  # where the hiptl files are saved

# input files
ptlfile = hp.File(PTLPATH+"snap_%03d.%d.hdf5" %(SNAPSHOT, CHUNK), 'r')

# output files
if IN_RS_SPACE:
    w = hp.File(OUTPATH+'ptlrs%d_%03d.%d.hdf5' %(BOX, SNAPSHOT, CHUNK), 'w')
else:
    w = hp.File(OUTPATH+'ptl%d_%03d.%d.hdf5' %(BOX, SNAPSHOT, CHUNK), 'w')

# getting author-defined constants (these COULD change but are not expected to)
GRID = (RES,RES,RES)

# getting the needed simulation-defined constants
head = dict(ptlfile['Header'].attrs)
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
BOXSIZE = head['BoxSize']/1e3 #Mpc/h
REDSHIFT = head['Redshift']
SCALE_FACTOR = head['Time']
DMPTL = head['MassTable'][1] *1e10/LITTLE_H # solar masses
ptltype = [0,1,4,5]
ptltype = ['PartType'+str(i) for i in ptltype]
nptl = head['NumPart_ThisFile']


# getting data
field = np.zeros(GRID, dtype=np.float32)
for p in ptltype:
    if '1' in p:
        mass = np.ones(nptl[1]) * DMPTL
    else:
        mass = ptlfile[p]['Masses'][:]*1e10/LITTLE_H # solar masses
    pos = ptlfile[p]['Coordinates'][:]/1e3 * SCALE_FACTOR # Mpc/h
    vel = ptlfile[p]['Velocities'][:] * np.sqrt(SCALE_FACTOR) # km/s

    # shifting the positions to redshift space
    if IN_RS_SPACE:
        pos = pos_redshift_space(pos, vel, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)

    # assigning them into the field using the Mass Assignment Scheme given
    CICW(pos,field,BOXSIZE,mass)

w.create_dataset("particles", data=field, compression="gzip", compression_opts=9)


w.close()
