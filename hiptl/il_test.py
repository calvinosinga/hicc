import numpy as np
import h5py as hp
import illustris_python as il
import time
from library_hicc.mas import CICW
from library_hicc.models import get_hiptl_models
from library_hicc.redshift_space import pos_redshift_space
import sys

TNG = '/lustre/cosinga/tng100/'
HIPATH = '/lustre/diemer/illustris/hih2/'  # where the hiptl files are saved

resos = (1024, 1200, 1400, 1600, 1800, 2048)

print("loading data")
fields = ['Masses','Coordinates','Velocities']
head = il.groupcat.loadHeader(TNG,99)
print("loading the header took %.3e")
data = il.snapshot.loadSubset(TNG,99,0,fields=fields)

nfiles = head['NumFiles']
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
BOXSIZE = head['BoxSize']/1e3 #Mpc/h
REDSHIFT = head['Redshift']
SCALE_FACTOR = head['Time']
w = hp.File("/lustre/cosinga/final_fields/il_test_neut_hydrogen.hdf5",'w')
nfrac = np.zeros_like(data['Density'])
cidx = 0
for i in range(nfiles):
    f = hp.File(HIPATH+'hih2_particles_099.%d.hdf5','r')
    temp = f['PartType0']['f_neutral_H'][:]
    nfrac[cidx:cidx+len(temp)]=temp
    f.close()
nfrac = np.where(nfrac >= 0, nfrac, np.zeros_like(nfrac))

mass = data['Masses'][:]*1e10/LITTLE_H*nfrac # solar masses
pos = data['Coordinates'][:]/1e3 # Mpc/h
vel = data['Velocities'][:] * np.sqrt(SCALE_FACTOR) # km/s
redpos = pos_redshift_space(pos, vel, BOXSIZE, LITTLE_H*100, REDSHIFT, 0)
for r in resos:
    field = np.zeros((r,r,r), dtype=np.float32)
    CICW(pos, field, BOXSIZE, mass)
    w.create_dataset("res=%d"%r, data=field, compression="gzip", compression_opts=9)
    field = np.zeros((r,r,r), dtype=np.float32)
    CICW(redpos, field, BOXSIZE, mass)
    w.create_dataset("rss_res=%d"%r, data=field, compression="gzip", compression_opts=9)

w.close()
