#!/usr/bin/env python3
"""
This file combines the given subfiles' fields for the ptl run to
create a final hdf5 file with all the needed information saved.
"""

import numpy as np
import h5py as hp
import sys

from numpy.core.defchararray import count

# setting author-defined variables (not expected to change)
BASE = '/lustre/cosinga/hiptl_output/'
FINAL = '/lustre/cosinga/final_fields/'
GRID = (2048,2048,2048)

# getting command-line input
PREFIX = sys.argv[1]
START = int(sys.argv[2])
END = int(sys.argv[3])
SNAPSHOT = int(sys.argv[4])
BOX = int(sys.argv[5])
STEP = int(sys.argv[6]) # tells if this is the 1st or 2nd step in combine process

# creating output file
w = hp.File(FINAL+'%s%d_%03d.final.hdf5'%(PREFIX, BOX, SNAPSHOT),'w')

# getting array of files to iterate over
if STEP == 0:
    w = hp.File(BASE+'%s%d_%03d.%d.%d.hdf5'%(PREFIX, BOX, SNAPSHOT, START, END),'w')
    filenos = np.arange(START, END)
    files = ['%s%d_%03d.%d.hdf5'%(PREFIX, BOX, SNAPSHOT, i) for i in filenos]
elif STEP == 1:
    w = hp.File(FINAL+'%s%d_%03d.final.hdf5'%(PREFIX, BOX, SNAPSHOT),'w')
    filenos = np.arange(START, END, 20)
    files = ['%s%d_%03d.%d.%d.hdf5'%(PREFIX, BOX, SNAPSHOT, i, i+20) for i in filenos]
else:
    raise ValueError("the STEP input must be 0 or 1")

print('first file: ' + files[0])
print('last file: ' + files[-1])

# array of files to iterate over
keylist = ((hp.File(BASE+files[0],'r')).keys()) # getting the keys from first file

total_counts = None
for k in keylist:
    total = np.zeros(GRID, dtype=np.float32)
    for i in files:
        f = hp.File(BASE+i,'r')
        if "count" in k: # since the number of bins should stay flexible
            if total_counts is None:
                total_counts = f[k][:]
            else:
                total_counts += f[k][:]
        else:
            total += f[k][:]
            print('found file %s for %s, adding to grid'%(i,k))
            print("new sum:" + str(np.sum(total)))
        f.close()
    w.create_dataset(k, data=total, compression="gzip", compression_opts=9)
w.create_dataset("bin_counts", data=total_counts)
w.close()
