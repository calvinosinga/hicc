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
from library_hicc.printer import Printer

# reading command line inputs
CHUNK = int(sys.argv[1])
SNAPSHOT = int(sys.argv[2])
BOX = int(sys.argv[3])
RES = int(sys.argv[4]) # resolution of the grid
AXIS = int(sys.argv[5])
GRID = (RES,RES,RES)

# defining needed paths
PTLPATH = '/lustre/cosinga/tng%d/snapdir_%03d/'%(BOX, SNAPSHOT) # where the ptl files are saved
OUTPATH = '/lustre/cosinga/HI-color/hicc/chunk_output/' # where to save the output
TREECOOL = '/lustre/cosinga/HI-color/TREECOOL_fg_dec11'
LOG = '/lustre/cosinga/HI-color/hicc/logs/v-n/'

# get input file
ptlfile = hp.File(PTLPATH+"snap_%03d.%d.hdf5" %(SNAPSHOT, CHUNK), 'r')

# output file
w = hp.File(OUTPATH+"v-n%d_%03d.%d.hdf5"%(BOX, SNAPSHOT, CHUNK),'w')
wrs = hp.File(OUTPATH+"v-nrs%d_%03d.%d.hdf5"%(BOX, SNAPSHOT, CHUNK),'w')
pnt = Printer(LOG+"v-n%d_%03d.%d.log"%(BOX, SNAPSHOT, CHUNK))

# getting the needed simulation-defined constants
head = dict(ptlfile['Header'].attrs)
SCALE_FACTOR = head['Time']
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
BOXSIZE = head['BoxSize']/1e3 * SCALE_FACTOR #Mpc/h
REDSHIFT = head['Redshift']

# creating HI using Paco's function
pnt.write("getting the HI from Paco's method...")
pos, MHI = vnhi(PTLPATH+"snap_%03d.%d.hdf5"%(SNAPSHOT, CHUNK),TREECOOL)

# returns positions in cMpc, mass in solar masses/h
MHI *= 1/LITTLE_H
pos *= SCALE_FACTOR # to convert from comoving

vel = ptlfile['PartType0']['Velocities'][:] * np.sqrt(SCALE_FACTOR) # km/s

pnt.write("creating the real-space grid...")
field = np.zeros(GRID, dtype=np.float32)
pnt.write("shifting the positions into redshift-space...")
rspos = pos_redshift_space(pos, vel, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)

pnt.write("placing the real-space data into the grid...")
CICW(pos,field,BOXSIZE,MHI)
pnt.write("saving the real-space data...")
w.create_dataset("v-n", data=field, compression="gzip", compression_opts=9)

pnt.write("creating new grid for redshift-space...")
field = np.zeros(GRID, dtype=np.float32)

pnt.write("placing redshift-space data into the grid...")
CICW(rspos,field,BOXSIZE,MHI)
pnt.write("saving the redshift-space data...")
wrs.create_dataset("v-n", data=field, compression="gzip", compression_opts=9)
wrs.close()
w.close()
ptlfile.close()
