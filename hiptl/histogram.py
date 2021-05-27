#!/usr/bin/env python3
import numpy as np
import h5py as hp
import sys

SNAPSHOT = int(sys.argv[1])
BOX = int(sys.argv[2])
NUMFILES = int(sys.argv[3])

PATH = "/lustre/cosinga/hiptl_output/"
f = hp.File(PATH + 'numdens%d_%03d.0.hdf5'%(BOX, SNAPSHOT),'r')
keylist = list(f.keys())
dims = f[keylist[0]][:].shape
w = hp.File('/lustre/cosinga/final_fields/numdens.hdf5','w')
for k in keylist:
    total = np.zeros(dims)
    for i in range(NUMFILES):
        f = hp.File(PATH + 'numdens%d_%03d.%d.hdf5'%(BOX, SNAPSHOT, i),'r')
        total += f[k][:]
    w.create_dataset(k, data=total)
w.close()

