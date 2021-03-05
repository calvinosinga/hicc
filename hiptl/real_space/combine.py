#!/bin/python3
"""
This file combines the given subfiles' fields for the hiptl run to
create a final hdf5 file with all the needed information saved.
"""

import numpy as np
import h5py as hp
import sys
from hicc_library.models import get_hiptl_models

# setting author-defined variables (not expected to change)
BASE = '/lustre/cosinga/hiptl_output/'
LOG = '/lustre/cosinga/hicc/hiptl/real_space/logs/'
models = get_hiptl_models()

# getting command-line input
START = sys.argv[1]
END = sys.argv[2]
SNAPSHOT = sys.argv[3]
BOX = sys.argv[4]

# opening files to write to
w = hp.File(BASE+'hiptl%s_%s.%s.%s.hdf5'%(BOX, SNAPSHOT, START, END),'w')
logfile = open(BASE+'combine%_%s.log'%(BOX, SNAPSHOT),'a')

# getting the filenames that we will be combining
filenos = np.arange(int(START), int(END))
files = ['hiptl%s_%s.%s.hdf5'%(BOX, SNAPSHOT, i) for i in filenos]

logfile.write('first file: ' + files[0]+'\n')
logfile.write('last file: ' + files[-1]+'\n')

for m in models:
    total = np.zeros((2048, 2048, 2048), dtype=np.float32)
    logfile.write("starting model "+m+'\n')
    logfile.write("current total sum %.4f"%(np.sum(total)))
    for i in files:
        f = hp.File(BASE+i,'r')
        total += f[m][:]
        logfile.write("new sum:" + str(np.sum(total))+"\n")
        f.close()
    w.create_dataset(m, data=total, compression="gzip", compression_opts=9)
w.close()
logfile.close()
