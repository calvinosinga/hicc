#!/usr/bin/env python3
"""
This file attempts to recreate Paco's fields
"""

# import statements
import numpy as np
import h5py as hp
import sys
from library_hicc.mas import CICW
from library_hicc.redshift_space import pos_redshift_space
from HI_library import HI_mass_from_Illustris_snap as vnhi

# reading command line inputs
CHUNK = int(sys.argv[1])
SNAPSHOT = int(sys.argv[2])
BOX = int(sys.argv[3])
RES = int(sys.argv[4]) # resolution of the grid
AXIS = int(sys.argv[5])
GRID = (RES,RES,RES)

# defining needed paths
PTLPATH = '/lustre/cosinga/tng%d/snapdir_%03d/'%(BOX, SNAPSHOT) # where the ptl files are saved
OUTPATH = '/lustre/cosinga/hiptl_output/' # where to save the output
TREECOOL = '/lustre/cosinga/fg20_treecool/TREECOOL_fg_dec11'

# get input file
ptlfile = hp.File(PTLPATH+"snap_%03d.%d.hdf5" %(SNAPSHOT, CHUNK), 'r')

# output file
w = hp.File(OUTPATH+"v-n_trial%d_%03d.%d.hdf5"%(BOX, SNAPSHOT, CHUNK))
wrs = hp.File(OUTPATH+"v-n_trialrs%d_%03d.%d.hdf5"%(BOX, SNAPSHOT, CHUNK))

# getting the needed simulation-defined constants
head = dict(ptlfile['Header'].attrs)
SCALE_FACTOR = head['Time']
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
BOXSIZE = head['BoxSize']/1e3 * SCALE_FACTOR #Mpc/h
REDSHIFT = head['Redshift']


# creating HI using Paco's function
pos, MHI = vnhi(PTLPATH+"snap_%03d.%d.hdf5"%(SNAPSHOT, CHUNK),TREECOOL)
# returns positions in cMpc, mass in solar masses/h
MHI *= 1/LITTLE_H
pos *= SCALE_FACTOR # to convert from comoving

vel = ptlfile['PartType0']['Velocities'][:] * np.sqrt(SCALE_FACTOR) # km/s

field = np.zeros(GRID, dtype=np.float32)
rspos = pos_redshift_space(pos, vel, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)

CICW(pos,field,BOXSIZE,MHI)
w.create_dataset("v-n", data=field, compression="gzip", compression_opts=9)

field = np.zeros(GRID, dtype=np.float32)
CICW(rspos,field,BOXSIZE,MHI)
wrs.create_dataset("v-n", data=field, compression="gzip", compression_opts=9)