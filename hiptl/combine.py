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

# this is giving memory error for some unknown reason, checking for why
import os, psutil
process = psutil.Process(os.getpid())

def mem_to_string():
    mem = process.memory_info().rss
    mem = mem/1e6
    return "%.3e"%mem
# getting command-line inputs
PREFIX = sys.argv[1]
START = int(sys.argv[2])
END = int(sys.argv[3])
SNAPSHOT = int(sys.argv[4])
BOX = int(sys.argv[5])
STEP = int(sys.argv[6]) # tells if this is the 1st or 2nd step in combine process

# setting author-defined variables (not expected to change)
models = get_hiptl_models()
BASE = "/lustre/cosinga/HI-color/chunk_output/"
FINAL = "/lustre/cosinga/HI-color/results/fields/snap_%03d/"%SNAPSHOT
LOG = '/lustre/cosinga/HI-color/hicc/logs/hiptl/'

# getting array of files to iterate over
if STEP == 0:
    w = hp.File(BASE+'%s%d_%03d.%d.%d.hdf5'%(PREFIX, BOX, SNAPSHOT, START, END),'w')
    pnt = Printer(LOG+'%s%d_%03d.%d.%d.log'%(PREFIX, BOX, SNAPSHOT, START, END))
    filenos = np.arange(START, END)
    files = ['%s%d_%03d.%d.hdf5'%(PREFIX, BOX, SNAPSHOT, i) for i in filenos]
elif STEP == 1:
    filenos = np.arange(START, END, 20)
    files = ['%s%d_%03d.%d.%d.hdf5'%(PREFIX, BOX, SNAPSHOT, i, i+20) for i in filenos]
else:
    raise ValueError("the STEP input must be 0 or 1")


f = hp.File(BASE+files[0], 'r')
keylist = list(f.keys())
GRID = f[models[0]][:].shape

if STEP==1:
    w = hp.File(FINAL+'%s%d_%03d.%dres.final.hdf5'%(PREFIX, BOX, SNAPSHOT,GRID[0]),'w')
    pnt = Printer(LOG+'%s%d_%03d.%dres.final.log'%(PREFIX, BOX, SNAPSHOT,GRID[0]))
    
pnt.write("just loaded input, about to start loop:"+mem_to_string())
# after input loading, has memory of 47.2MB
# sum each model's grid individually
for m in models:
    pnt.write("starting model "+m+"\tmem="+mem_to_string())
    total = np.zeros(GRID, dtype=np.float32)

    pnt.writeTab("created empty array:"+mem_to_string())
    # now using 3.5 GB
    for i in files:
        # it is expected that the last job will find nonexistent files
        try:
            f = hp.File(BASE+i,'r')
        except IOError:
            pnt.writeTab('did not find the file %s'%i)
        else:
            pnt.writeTab('found file %s, adding to grid'%i)
            pnt.writeTab("just loaded file:"+mem_to_string())
            total += f[m][:]
            pnt.writeTab("added to running total:"+mem_to_string())
            pnt.writeTab("new sum:" + str(np.sum(total)))
            f.close()
    w.create_dataset(m, data=total, compression="gzip", compression_opts=9)
    pnt.writeTab("saved total to output file"+mem_to_string())
w.close()
