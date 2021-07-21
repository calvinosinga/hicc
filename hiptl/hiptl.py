#!/usr/bin/env python3
"""
This file takes each of the 448 subfiles for the particle-based hydrogen prescriptions
and creates a 2048^3 grid of each of them.
"""
# import statements
import numpy as np
import h5py as hp
import sys
from library_hicc.mas import CICW
from library_hicc.models import get_hiptl_models
from library_hicc.redshift_space import pos_redshift_space
from library_hicc.printer import Printer

# reading command line inputs
CHUNK = int(sys.argv[1])
SNAPSHOT = int(sys.argv[2])
BOX = int(sys.argv[3])
RES = int(sys.argv[4]) # resolution of the grid
AXIS = int(sys.argv[5]) # if -1, not in redshift space
GRID = (RES,RES,RES)

# defining needed paths
HIPATH = '/lustre/diemer/illustris/hih2/'  # where the hiptl files are saved
PTLPATH = '/lustre/cosinga/tng%d/snapdir_%03d/'%(BOX, SNAPSHOT) # where the ptl files are saved
OUTPATH = '/lustre/cosinga/chunk_output/' # where to save the output
LOG = '/lustre/cosinga/HI-color/hicc/logs/hiptl/'

# input files
hih2file = hp.File(HIPATH+"hih2_particles_%03d.%d.hdf5" %(SNAPSHOT, CHUNK), 'r')
ptlfile = hp.File(PTLPATH+"snap_%03d.%d.hdf5" %(SNAPSHOT, CHUNK), 'r')

# output files
wrs = hp.File(OUTPATH+'hiptlrs%d_%03d.%d.hdf5' %(BOX, SNAPSHOT, CHUNK), 'w')
w = hp.File(OUTPATH+'hiptl%d_%03d.%d.hdf5' %(BOX, SNAPSHOT, CHUNK), 'w')
pnt = Printer(LOG+'hiptl%d_%03d.%d.log' %(BOX, SNAPSHOT, CHUNK))

# getting author-defined constants (these COULD change but are not expected to)
models = get_hiptl_models()

# getting the needed simulation-defined constants
head = dict(ptlfile['Header'].attrs)
SCALE = head['Time'] # scale factor
BOXSIZE = head['BoxSize']/1e3 * SCALE #Mpc/h
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
REDSHIFT = head['Redshift']

# getting data
mass = ptlfile['PartType0']['Masses'][:]*1e10/LITTLE_H # solar masses
pos = ptlfile['PartType0']['Coordinates'][:]/1e3 * SCALE # Mpc/h
vel = ptlfile['PartType0']['Velocities'][:] * np.sqrt(SCALE) # km/s
f_neut_h = hih2file['PartType0']['f_neutral_H'][:]

for m in models:
    pnt.write('starting to calculate grid for %s'%m)
    field = np.zeros(GRID, dtype=np.float32)

    # getting the HI mass data
    h2_frac = hih2file['PartType0']['f_mol_'+m][:]
    masshi = (1-h2_frac)*f_neut_h*mass
    pnt.writeTab('finished calculating the hi mass, has sum %.3e'%np.sum(masshi))

    # f_neut_h is negative where the model isn't defined, removing negative masses.
    masshi = np.where(masshi >= 0, masshi, np.zeros(masshi.shape, dtype=np.float32))
    masshi = masshi.astype('float32')
    pnt.writeTab('finished removing the negative fractions, now has sum %.3e'%np.sum(masshi))
    
    # assigning them into the field using the Mass Assignment Scheme given
    pnt.writeTab('now assigning the real-space HI to a grid')
    CICW(pos,field,BOXSIZE,masshi)
    pnt.writeTab('saving the grid...')
    w.create_dataset(m, data=field, compression="gzip", compression_opts=9)

    # now doing the same for redshift space
    pnt.writeTab('creating a new grid for redshift-space')
    field = np.zeros(GRID, dtype=np.float32)

    # shifting the positions to redshift space   
    pnt.writeTab('shifting positions to redshift space ...')
    rspos = pos_redshift_space(pos, vel, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)

    pnt.writeTab('now assigning the redshift-space HI to a grid')
    CICW(rspos, field, BOXSIZE, masshi)
    pnt.writeTab('saving the redshift-space grid')
    wrs.create_dataset(m, data=field, compression="gzip", compression_opts=9)

w.close()
wrs.close()