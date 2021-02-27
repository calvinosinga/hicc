#!/usr/bin/env python3
"""
This file takes each of the 448 subfiles for the particle-based hydrogen prescriptions
and creates a 2048^3 grid of each of them.

In addition, it'll save (for the purpose of visualization) the exact positions of every
particle within some distance of a subhaloes of interest (SOI), like a Milky Way-like 
galaxy, a red galaxy and a blue galaxy.
"""
# import statements
import numpy as np
import h5py as hp
import sys
import MAS_library as masl

# reading command line inputs
CHUNK = sys.argv[1]
SNAPSHOT = sys.argv[2]

# defining needed paths
HIPATH = '/lustre/diemer/illustris/hih2/'  # where the hiptl files are saved
PTLPATH = '/lustre/cosinga/ptl%s/'%SNAPSHOT # where the ptl files are saved
OUTPATH = '/lustre/cosinga/hiptl_output/' # where to save the output
SUBPATH = '/lustre/cosinga/subhalo%s/'%SNAPSHOT # to find positions of SOI

# input files
hih2file = hp.File(HIPATH+"hih2_particles_0%s.%s.hdf5" %(SNAPSHOT, CHUNK), 'r')
ptlfile = hp.File(PTLPATH+"snap_0%s.%s.hdf5" %(SNAPSHOT, CHUNK), 'r')

# output files
w = hp.File(OUTPATH+'hiptl_%s.%s.hdf5' %(SNAPSHOT, CHUNK), 'w')

# getting author-defined constants (these COULD change but are not expected to)
MAS = 'CIC'
GRID = (2048,2048,2048)

# getting the needed simulation-defined constants
head = dict(ptlfile['Header'].attrs)
LITTLE_H = head['HubbleParam']
BOXSIZE = head['BoxSize']/1e3 #Mpc/h


# getting data
mass = ptlfile['PartType0']['Masses'][:]*1e10/LITTLE_H #~50 MB, in solar masses
pos = ptlfile['PartType0']['CenterOfMass'][:]/1e3 # ~150 MB, Mpc/h
f_neut_h = hih2file['PartType0']['f_neutral_H'][:] #~100 MB

# loop over each model
models = ['GD14', 'GK11', 'K13', 'S14']
for m in models:
    field = np.zeros(GRID, dtype=np.float32)

    # getting the HI mass data
    h2_frac = hih2file['PartType0']['f_mol_'+m][:] # ~100 MB
    masshi = (1-h2_frac)*f_neut_h*mass # ~100 MB

    # f_neut_h is negative where the model isn't defined, removing negative masses.
    masshi = np.where(masshi >= 0, masshi, np.zeros(masshi.shape, dtype=np.float32))
    masshi = masshi.astype('float32')

    # assigning them into the field using the Mass Assignment Scheme given
    masl.MA(pos,field,BOXSIZE,MAS,masshi)
    w.create_dataset(m, data=field, compression="gzip", compression_opts=9)
    
    #TODO: save the nearby HI here
w.close()
