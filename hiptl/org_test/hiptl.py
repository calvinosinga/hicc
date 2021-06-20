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
# import redshift_space_library as rsl
from library_hicc.redshift_space import pos_redshift_space

# reading command line inputs
START = int(sys.argv[1])
END = int(sys.argv[2])
SNAPSHOT = int(sys.argv[3])
BOX = int(sys.argv[4])
RES = int(sys.argv[5]) # resolution of the grid
AXIS = int(sys.argv[6]) # if -1, not in redshift space
GRID = (RES,RES,RES)

# defining needed paths
HIPATH = '/lustre/diemer/illustris/hih2/'  # where the hiptl files are saved
PTLPATH = '/lustre/cosinga/tng%d/snapdir_%03d/'%(BOX, SNAPSHOT) # where the ptl files are saved
OUTPATH = '/lustre/cosinga/hiptl_output/' # where to save the output

# output files
w = hp.File(OUTPATH+'hiptlrs%d_%03d.%d.%d.hdf5' %(BOX, SNAPSHOT, START,END), 'w')
wrs = hp.File(OUTPATH+'hiptl%d_%03d.%d.%d.hdf5' %(BOX, SNAPSHOT, START, END), 'w')

# files to iterate over
filenos = np.arange(START, END)

# files to access
hih2files = ["%shih2_particles_%03d.%d.hdf5" %(HIPATH, SNAPSHOT, i) for i in filenos]
ptlfiles = [PTLPATH+"snap_%03d.%d.hdf5" %(SNAPSHOT, i) for i in filenos]

# getting the needed simulation-defined constants
tempptl = hp.File(ptlfiles[0], 'r')
head = dict(tempptl['Header'].attrs)
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
BOXSIZE = head['BoxSize']/1e3 #Mpc/h
REDSHIFT = head['Redshift']
SCALE_FACTOR = head['Time']
del tempptl

# models to go over
models = get_hiptl_models()

# particle shorthand
p = 'PartType0'
for m in models:
    # we iterate over the models first to save on memory - it means we have to loop
    # over the files several times but 16 GB per field x 4 models x 2 spaces gives more 
    # memory than one node has so we need to sacrifice run time
    print('starting to calculate grids for %s'%m)
    total = np.zeros(GRID, dtype=np.float32)
    totalrs = np.zeros(GRID, dtype=np.float32)
    for i in range(len(filenos)):
        # handle IOerrors when trying to open files (expected that some files 
        # will not exist in the last job)
        try:
            hih2f = hp.File(hih2files[i],'r')
            ptlf = hp.File(ptlfiles[i],'r')
        except IOError:
            print('did not find the file %s,%s'%(hih2files[i],ptlfiles[i]))
        else:
            # getting data
            mass = ptlf[p]['Masses'][:]*1e10/LITTLE_H # solar masses
            pos = ptlf[p]['Coordinates'][:]/1e3 * SCALE_FACTOR # Mpc/h
            vel = ptlf[p]['Velocities'][:] * np.sqrt(SCALE_FACTOR) # km/s
            f_neut_h = hih2f[p]['f_neutral_H'][:]
            
            h2_frac = hih2f[p]['f_mol_'+m][:]
            masshi = (1-h2_frac)*f_neut_h*h2_frac

            # removing negative masses (-1 is given where hih2 frac is undefined)
            masshi = np.where(masshi >= 0, masshi, np.zeros(masshi.shape, dtype=np.float32))
            masshi = masshi.astype('float32')
            print('finished removing negative masses, now has sum %.3e'%np.sum(masshi))

            # assigning real-space to the field using CICW mass-assignment scheme
            CICW(pos,total, BOXSIZE, masshi)

            # now putting the positions into redshift-space
            pos = pos_redshift_space(pos, vel, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)
            CICW(pos, totalrs, BOXSIZE, masshi)
            hih2f.close()
            ptlf.close()
    w.create_dataset(m, data=total, compression="gzip", compression_opts=9)
    wrs.create_dataset(m, data=totalrs, compression="gzip", compression_opts=9)
w.close()
wrs.close()
        




            
