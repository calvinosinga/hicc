#!/usr/bin/env python3
"""
This file combines the given subfiles' fields for the ptl run to
create a final hdf5 file with all the needed information saved.
"""

import numpy as np
import h5py as hp
import sys

# setting author-defined variables (not expected to change)
BASE = '/lustre/cosinga/hiptl_output/'
FINAL = '/lustre/cosinga/final_fields/'
models = get_hiptl_models()

# getting command-line input
PREFIX = sys.argv[1]
SNAPSHOT = int(sys.argv[2])
BOX = int(sys.argv[3])
FILENUM = int(sys.argv[4]) # tells how many files there are to combine

# creating output file
w = hp.File(FINAL+'%s%d_%03d.final.hdf5'%(PREFIX, BOX, SNAPSHOT),'w')

# array of files to iterate over

filenos = np.arange(0,FILENUM)
files = ['%s%d_%03d.%d.hdf5'%(PREFIX, BOX, SNAPSHOT, i) for i in filenos]
print('first file: ' + files[0])
print('last file: ' + files[-1])

total = np.zeros((2048, 2048, 2048), dtype=np.float32)

for i in files:
    f = hp.File(BASE+i,'r')
    total += f[m][:]
    print('found file %s, adding to grid'%i)
    print("new sum:" + str(np.sum(total)))
    f.close()
w.create_dataset("particles", data=total, compression="gzip", compression_opts=9)
w.close()
