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
if len(sys.argv) > 4:
    AXIS = int(sys.argv[4])
    IN_RS_SPACE = True
else:
    IN_RS_SPACE = False

# defining needed paths
PTLPATH = '/lustre/cosinga/tng%d/snapdir_%03d/'%(BOX, SNAPSHOT) # where the ptl files are saved
OUTPATH = '/lustre/cosinga/hiptl_output/' # where to save the output
HIPATH = '/lustre/diemer/illustris/hih2/'  # where the hiptl files are saved

# input files
ptlfile = hp.File(PTLPATH+"snap_%03d.%d.hdf5" %(SNAPSHOT, CHUNK), 'r')
hih2file = hp.File(HIPATH+"hih2_particles_%03d.%d.hdf5" %(SNAPSHOT, CHUNK), 'r')

# output files
if IN_RS_SPACE:
    w = hp.File(OUTPATH+'ptlrs%d_%03d.%d.hdf5' %(BOX, SNAPSHOT, CHUNK), 'w')
else:
    w = hp.File(OUTPATH+'ptl%d_%03d.%d.hdf5' %(BOX, SNAPSHOT, CHUNK), 'w')

# getting author-defined constants (these COULD change but are not expected to)
GRID = (2048,2048,2048)

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


# # getting data
# field = np.zeros(GRID, dtype=np.float32)
# for p in ptltype:
    

#     if '1' in p:
#         mass = np.ones(nptl[1]) * DMPTL
#     else:
#         mass = ptlfile[p]['Masses'][:]*1e10/LITTLE_H # solar masses
#     pos = ptlfile[p]['Coordinates'][:]/1e3 # Mpc/h
#     vel = ptlfile[p]['Velocities'][:] * np.sqrt(SCALE_FACTOR) # km/s

#     # shifting the positions to redshift space
#     if IN_RS_SPACE:
#         pos = pos_redshift_space(pos, vel, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)

#     # assigning them into the field using the Mass Assignment Scheme given
#     CICW(pos,field,BOXSIZE,mass)

# w.create_dataset("particles", data=field, compression="gzip", compression_opts=9)


# testing the different hydrogen densities

dendec =[0,1e-4,1e-2,1,100,np.inf]
counts = np.zeros_like(dendec)
p = 'PartType0'
mass = ptlfile[p]['Masses'][:]*1e10/LITTLE_H # solar masses
density = ptlfile[p]['Density'][:]*1e10/LITTLE_H # solar masses per kpc/h cubed
# tothyd = ptlfile[p]['GFM_Metals'][:,0] # fraction of cell that is hydrogen
f_neut_h = hih2file['PartType0']['f_neutral_H'][:]
# replacing the cells with negative fractions with zero
f_neut_h = np.where(f_neut_h>=0, f_neut_h, np.zeros_like(f_neut_h))
kpctocm = 3.086e21
smtog = 1.989e33
m_p=1.673e-24
fac = 1/kpctocm**3*LITTLE_H**3*smtog/m_p

n_h = density*f_neut_h*fac
pos = ptlfile[p]['Coordinates'][:]/1e3 # Mpc/h
vel = ptlfile[p]['Velocities'][:] * np.sqrt(SCALE_FACTOR) # km/s
print("now binning according to number density, there are " + str(n_h.shape)+ " cells.")
for d in range(len(dendec)-1):
    field = np.zeros(GRID, dtype=np.float32)
    lo = dendec[d]
    hi = dendec[d+1]
    print("getting low mask: %.3e"%lo)
    mask1 = n_h >= lo
    print("low mask has %d cells in it"%np.sum(mask1))
    print("getting high mask: %.3e"%hi)
    mask2 = n_h < hi
    print("high mask has %d cells in it"%np.sum(mask2))
    mask = mask1 & mask2
    counts[d] = np.sum(mask)
    print("in the <%.3e bin there are %d cells, has average %.4e"%(hi, counts[d], np.mean(n_h[mask])))
    posmask = pos[mask,:]
    massmask = mass[mask]
    if IN_RS_SPACE:
        velmask = vel[mask,:]
        posmask = pos_redshift_space(posmask, velmask, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)
    CICW(posmask,field,BOXSIZE,massmask)
    w.create_dataset("<%.3e"%dendec[d], data=field, compression="gzip", compression_opts=9)
w.create_dataset("bin_counts",data=counts)
w.close()
