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
BASE = '/lustre/cosinga/hiptl_output/'
FINAL = '/lustre/cosinga/final_fields/'
models = get_hiptl_models()

# getting command-line input
START = int(sys.argv[1])
END = int(sys.argv[2])
SNAPSHOT = int(sys.argv[3])
BOX = int(sys.argv[4])
STEP = int(sys.argv[5]) # tells if this is the 1st or 2nd step in combine process

# opening files to write to, getting the filenames of the files we are combining
if STEP == 0:
    w = hp.File(BASE+'hiptl%d_%03d.%d.%d.hdf5'%(BOX, SNAPSHOT, START, END),'w')
    filenos = np.arange(START, END)
    files = ['hiptl%d_%03d.%d.hdf5'%(BOX, SNAPSHOT, i) for i in filenos]
elif STEP == 1:
    w = hp.File(FINAL+'hiptl%d_%03d.final.hdf5'%(BOX, SNAPSHOT),'w')
    filenos = np.arange(START, END, 20)
    files = ['hiptl%d_%03d.%d.%d.hdf5'%(BOX, SNAPSHOT, i, i+20) for i in filenos]
else:
    raise ValueError("the STEP input must be 0 or 1")

print('first file: ' + files[0]+'\n')
print('last file: ' + files[-1]+'\n')

# sum each model's grid individually
for m in models:
    total = np.zeros((2048, 2048, 2048), dtype=np.float32)
    print("starting model "+m+'\n')
    print("current total sum %.4f"%(np.sum(total)))
    for i in files:
        # it is expected that the last job will have nonexistant files
        try:
            f = hp.File(BASE+i,'r')
        except IOError:
            print('did not find the file %s\n'%i)
        else:
            total += f[m][:]
            print("new sum:" + str(np.sum(total))+"\n")
            f.close()
    w.create_dataset(m, data=total, compression="gzip", compression_opts=9)
w.close()
