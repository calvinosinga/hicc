#!/usr/bin/env python3
"""
This file combines the given subfiles' fields for the hiptl run to
create a final hdf5 file with all the needed information saved.
"""

import numpy as np
import h5py as hp
import sys
from library_hicc.models import get_hiptl_models
from library_hicc.printer import Printer

# getting command-line inputs
# for whatever reason this is giving an error
print(sys.argv)
PREFIX = sys.argv[1]
START = int(sys.argv[2])
END = int(sys.argv[3])
SNAPSHOT = int(sys.argv[4])
BOX = int(sys.argv[5])
STEP = int(sys.argv[6]) # tells if this is the 1st or 2nd step in combine process

# setting author-defined variables (not expected to change)
BASE = "/lustre/cosinga/HI-color/chunk_output/"
FINAL = "/lustre/cosinga/HI-color/results/fields/snap_%03d/"%SNAPSHOT
LOG = '/lustre/cosinga/HI-color/hicc/logs/v-n/'

# getting array of files to iterate over
if STEP == 0:
    w = hp.File(BASE+'%s%d_%03d.%d.%d.hdf5'%(PREFIX, BOX, SNAPSHOT, START, END),'w')
    pnt = Printer(LOG+'%s%d_%03d.%d.%d.log'%(PREFIX, BOX, SNAPSHOT, START, END))
    filenos = np.arange(START, END)
    files = ['%s%d_%03d.%d.hdf5'%(PREFIX, BOX, SNAPSHOT, i) for i in filenos]
elif STEP == 1:
    w = hp.File(FINAL+'%s%d_%03d.final.hdf5'%(PREFIX, BOX, SNAPSHOT),'w')
    pnt = Printer(LOG+'%s%d_%03d.final.log'%(PREFIX, BOX, SNAPSHOT))
    filenos = np.arange(START, END, 20)
    files = ['%s%d_%03d.%d.%d.hdf5'%(PREFIX, BOX, SNAPSHOT, i, i+20) for i in filenos]
else:
    raise ValueError("the STEP input must be 0 or 1")

pnt.write('first file: ' + files[0])
pnt.write('last file: ' + files[-1])
f = hp.File(BASE+files[0], 'r')
keylist = list(f.keys())
GRID = f[keylist[0]][:].shape
# sum each model's grid individually
for key in keylist:
    total = np.zeros(GRID, dtype=np.float32)
    pnt.write("starting model "+key)
    for i in files:
        # it is expected that the last job will find nonexistent files
        try:
            f = hp.File(BASE+i,'r')
        except IOError:
            pnt.writeTab('did not find the file %s'%i)
        else:
            total += f[key][:]
            pnt.writeTab('found file %s, adding to grid'%i)
            pnt.writeTab("new sum:" + str(np.sum(total)))
            f.close()
    w.create_dataset(key, data=total, compression="gzip", compression_opts=9)
w.close()
