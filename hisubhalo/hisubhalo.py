import numpy as np
import h5py as hp
import sys
import MAS_library as masl

# getting constants
grid = (2048,2048,2048)
SNAPSHOT = sys.argv[1]
BOXSIZE = 75.0 #Mpc/h
HOME = '/lustre/cosinga/subhalo'+str(SNAPSHOT)+'/'
SAVE = '/lustre/cosinga/subhalo_output/'
MAS = sys.argv[2]
###################################
# getting data and opening a log file
logfile = open(SAVE+'hisubhalo_log'+str(SNAPSHOT)+'.txt', 'a')
w = hp.File(SAVE+'hisubhalo_'+str(SNAPSHOT)+'.final.hdf5', 'w')
f = hp.File(HOME+'hih2_galaxy_0'+str(SNAPSHOT)+'.hdf5','r')
idfile = hp.File(HOME+"id_pos"+str(SNAPSHOT)+".hdf5",'r')
pos = idfile['coordinates'][:]/1e3 # Mpc/h
keys = list(f.keys())
models = []
for k in keys:
    if 'm_hi' in k:
        models.append(k)
logfile.write("the models used are: "+str(models)+'\n')

# loop over 9 models for HI
for m in models:
    field = np.zeros(grid, dtype=np.float32)
    mass = f[m][:] # already in solar masses
    mass = mass.astype(np.float32)
    pos = pos.astype(np.float32)
    masl.MA(pos, field, BOXSIZE, MAS, mass)
    w.create_dataset(m, data=field, compression="gzip", compression_opts=9)

w.close()
