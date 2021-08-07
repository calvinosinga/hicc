#!/usr/bin/env python3
"""
This file takes each of the 448 subfiles for the particle subfiles and creates a 2048^3
grid for them.
"""
# import statements
from functools import wraps
import numpy as np
import h5py as hp
import sys
from library_hicc.mas import CICW
from library_hicc.printer import Printer
from library_hicc.redshift_space import pos_redshift_space
import library_hicc.color as coldef
# reading command line inputs
CHUNK = int(sys.argv[1])
SNAPSHOT = int(sys.argv[2])
BOX = int(sys.argv[3])
RES = int(sys.argv[4]) # resolution of the grid
AXIS = int(sys.argv[5]) # if -1, not in redshift space

# defining needed paths
PTLPATH = '/lustre/cosinga/tng%d/snapdir_%03d/'%(BOX, SNAPSHOT) # where the ptl files are saved
OUTPATH = '/lustre/cosinga/HI-color/chunk_output/' # where to save the output
HIPATH = '/lustre/diemer/illustris/hih2/'  # where the hiptl files are saved
LOG = '/lustre/cosinga/HI-color/hicc/logs/ptl/'
GCPATH = '/lustre/cosinga/tng%d/groups_%03d/'%(BOX,SNAPSHOT)

# input files
ptlfile = hp.File(PTLPATH+"snap_%03d.%d.hdf5" %(SNAPSHOT, CHUNK), 'r')

# output files
wrs = hp.File(OUTPATH+'ptlrs%d_%03d.%d.hdf5' %(BOX, SNAPSHOT, CHUNK), 'w')
w = hp.File(OUTPATH+'ptl%d_%03d.%d.hdf5' %(BOX, SNAPSHOT, CHUNK), 'w')
pnt = Printer(LOG+'ptl%d_%03d.%d.log' %(BOX, SNAPSHOT, CHUNK))

# getting author-defined constants (these COULD change but are not expected to)
GRID = (RES,RES,RES)

# getting the needed simulation-defined constants
pnt.write("loading sim info...")
head = dict(ptlfile['Header'].attrs)
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
SCALE = head['Time'] # scale factor
BOXSIZE = head['BoxSize']/1e3 * SCALE #Mpc/h
REDSHIFT = head['Redshift']

DMPTL = head['MassTable'][1] *1e10/LITTLE_H # solar masses
ptltype = [0,1,4,5]
ptltype = ['PartType'+str(i) for i in ptltype]
nptl = head['NumPart_ThisFile']


# getting data
pnt.write('creating grid...')
field = np.zeros(GRID, dtype=np.float32)

pnt.write('starting loop for real-space...')
for p in ptltype:
    pnt.writeTab("getting data for %s"%p)
    if '1' in p:
        mass = np.ones(nptl[1], dtype=np.float32) * DMPTL
    else:
        mass = ptlfile[p]['Masses'][:]*1e10/LITTLE_H # solar masses
    pos = ptlfile[p]['Coordinates'][:]/1e3 * SCALE # Mpc/h

    # assigning them into the field using the Mass Assignment Scheme given
    pnt.writeTab("assigning to real-space grid...")
    CICW(pos,field,BOXSIZE,mass)

pnt.write("saving real-space grid...")
w.create_dataset("particles", data=field, compression="gzip", compression_opts=9)

# creating grids for galaxy-ptl
names = ['blue_ptl', 'red_ptl','unresolved_ptl']
subh = hp.File(GCPATH+'fof_subhalo_tab_%03d.%d.hdf5'%(SNAPSHOT,CHUNK),'r')
try:
    subhpos = subh['SubhaloPos'][:]
    subhphoto = subh['SubhaloStellarPhotometrics'][:]
    subhmass = subh['SubhaloMassType'][:]
    found_file = True
except KeyError:
    pnt.write("did not find any information in the file for %d"%CHUNK)
    found_file=False
field = np.zeros(GRID, dtype=np.float32)
if found_file:
    
w.close()

pnt.write('creating new grid for redshift space...')
field = np.zeros(GRID, dtype=np.float32)

pnt.write('starting loop for redshift-space...')
for p in ptltype:
    pnt.writeTab("getting data for %s"%p)
    if '1' in p:
        mass = np.ones(nptl[1]) * DMPTL
    else:
        mass = ptlfile[p]['Masses'][:]*1e10/LITTLE_H # solar masses
    pos = ptlfile[p]['Coordinates'][:]/1e3 * SCALE # Mpc/h
    vel = ptlfile[p]['Velocities'][:] * np.sqrt(SCALE) # km/s
    pnt.writeTab("moving positions to redshift space...")
    pos = pos_redshift_space(pos, vel, BOXSIZE, LITTLE_H*100, REDSHIFT, AXIS)
    # assigning them into the field using the Mass Assignment Scheme given
    pnt.writeTab("assigning to redshift-space grid...")
    CICW(pos,field,BOXSIZE,mass)

pnt.write("saving redshift-space grid...")
wrs.create_dataset("particles", data=field, compression="gzip", compression_opts=9)
wrs.close()