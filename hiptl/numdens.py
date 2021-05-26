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

w = hp.File(OUTPATH+'numdens%d_%03d.%d.hdf5' %(BOX, SNAPSHOT, CHUNK), 'w')

# getting author-defined constants (these COULD change but are not expected to)
GRID = (2048,2048,2048)
models = get_hiptl_models()

# getting the needed simulation-defined constants
head = dict(ptlfile['Header'].attrs)
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
BOXSIZE = head['BoxSize']/1e3 #Mpc/h
REDSHIFT = head['Redshift']
SCALE_FACTOR = head['Time']

p = "PartType0"
density = ptlfile[p]['Density'][:]*1e10/LITTLE_H # solar masses per kpc/h cubed
tothyd = ptlfile[p]['GFM_Metals'][:,0] # fraction of cell that is hydrogen
nha = ptlfile[p]["NeutralHydrogenAbundance"][:] # fraction of hydrogen that is neutral
kpctocm = 3.086e21
smtog = 1.989e33
m_p=1.673e-24

f_neut_H = hih2file[p]['f_neutral_H'][:] # fraction of total cell that is neutral hydrogen

# replacing the cells with negative fractions with zero
f_neut_H = np.where(f_neut_H>=0, f_neut_H, np.zeros_like(f_neut_H))

bins = [0,1e-4,1e-3,1e-2,1e-1,1,10,100,np.inf]

# converting density to g/cm^3
density = 1 / kpctocm**3 * LITTLE_H**3 * smtog * density

# using simple density to number density conversion rho = 1.4 mp nH

nh = density / m_p / 1.4
hist = np.histogram(nh, bins=bins)
w.create_dataset("simple", data=hist[0])

# using snapshot's definition of hydrogen abundance
nh = density * tothyd / m_p
hist = np.histogram(nh, bins=bins)
w.create_dataset("all_hydrogen", data=hist[0])

# using the snapshot's definition of nha without sfr adjustment
nh = density * tothyd * nha / m_p
hist = np.histogram(nh, bins=bins)
w.create_dataset("neutral_hydrogen_with_sfc", data=hist[0])

# using the Benedikt's definition of neutral hydrogen
nh = density * f_neut_H / m_p
hist = np.histogram(nh, bins=bins)
w.create_dataset("hih2_neutral_hydrogen", data=hist[0])

w.close()
