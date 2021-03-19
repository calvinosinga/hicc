#!/usr/bin/env python3
"""
This file takes each of the 448 subfiles for the particle-based hydrogen prescriptions
and creates a 2048^3 grid of each of them.
"""
# import statements
import numpy as np
import h5py as hp
import sys
# import MAS_library as masl
from library_hicc.mas import CICW
from library_hicc.models import get_hiptl_models
import redshift_space_library as rsl

# reading command line inputs
CHUNK = int(sys.argv[1])
SNAPSHOT = int(sys.argv[2])
BOX = int(sys.argv[3])
if len(sys.argv) > 4:
    AXIS = int(sys.argv[4])
    IN_RS_SPACE = True
else:
    IN_RS_SPACE = False

# defining needed paths
HIPATH = '/lustre/diemer/illustris/hih2/'  # where the hiptl files are saved
PTLPATH = '/lustre/cosinga/tng%d/snapdir_%03d/'%(BOX, SNAPSHOT) # where the ptl files are saved
OUTPATH = '/lustre/cosinga/hiptl_output/' # where to save the output

# input files
hih2file = hp.File(HIPATH+"hih2_particles_%03d.%d.hdf5" %(SNAPSHOT, CHUNK), 'r')
ptlfile = hp.File(PTLPATH+"snap_%03d.%d.hdf5" %(SNAPSHOT, CHUNK), 'r')

# output files
w = hp.File(OUTPATH+'hiptlrs%d_%03d.%d.hdf5' %(BOX, SNAPSHOT, CHUNK), 'w')

# getting author-defined constants (these COULD change but are not expected to)
MAS = 'CIC'
GRID = (2048,2048,2048)
models = get_hiptl_models()

# getting the needed simulation-defined constants
head = dict(ptlfile['Header'].attrs)
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
BOXSIZE = head['BoxSize']/1e3 #Mpc/h
REDSHIFT = head['Redshift']
SCALE_FACTOR = head['Time']

# getting data
mass = ptlfile['PartType0']['Masses'][:]*1e10/LITTLE_H # solar masses
pos = ptlfile['PartType0']['CenterOfMass'][:]/1e3 # Mpc/h
vel = ptlfile['PartType0']['Velocities'][:] * np.sqrt(SCALE_FACTOR) # km/s
f_neut_h = hih2file['PartType0']['f_neutral_H'][:]

for m in models:
    print('starting to calculate grid for %s'%m)
    field = np.zeros(GRID, dtype=np.float32)
    # shifting the positions to redshift space
    if IN_RS_SPACE:
        rsl.pos_redshift_space(pos, vel, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)

    # getting the HI mass data
    h2_frac = hih2file['PartType0']['f_mol_'+m][:]
    masshi = (1-h2_frac)*f_neut_h*mass
    print('finished calculating the hi mass, has sum %.3e'%np.sum(masshi))

    # f_neut_h is negative where the model isn't defined, removing negative masses.
    masshi = np.where(masshi >= 0, masshi, np.zeros(masshi.shape, dtype=np.float32))
    masshi = masshi.astype('float32')
    print('finished removing the negative fractions, now has sum %.3e'%np.sum(masshi))
    
    # assigning them into the field using the Mass Assignment Scheme given
    CICW(pos,field,BOXSIZE,masshi)
    w.create_dataset(m, data=field, compression="gzip", compression_opts=9)

w.close()

# checking if the output file was formatted as expected