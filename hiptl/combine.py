#!/usr/bin/env python3
"""
This file combines the given subfiles' fields for the hiptl run to
create a final hdf5 file with all the needed information saved.
"""

import numpy as np
import h5py as hp
import sys
from library_hicc.models import get_hiptl_models

# setting author-defined variables (not expected to change)
models = get_hiptl_models()

# getting command-line inputs
PREFIX = sys.argv[1]
START = int(sys.argv[2])
END = int(sys.argv[3])
SNAPSHOT = int(sys.argv[4])
BOX = int(sys.argv[5])
STEP = int(sys.argv[6]) # tells if this is the 1st or 2nd step in combine process

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

# creating output file
w = hp.File(FINAL+'%s%d_%03d.final.hdf5'%(PREFIX, BOX, SNAPSHOT),'w')

# sum each model's grid individually
for m in models:
    total = np.zeros((2048, 2048, 2048), dtype=np.float32)
    print("starting model "+m)
    for i in files:
        # it is expected that the last job will find nonexistent files
        try:
            f = hp.File(BASE+i,'r')
        except IOError:
            print('did not find the file %s'%i)
        else:
            total += f[m][:]
            print('found file %s, adding to grid'%i)
            print("new sum:" + str(np.sum(total)))
            f.close()
    w.create_dataset(m, data=total, compression="gzip", compression_opts=9)
w.close()
