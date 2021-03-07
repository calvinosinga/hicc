#!
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
LOG = '/lustre/cosinga/hicc/hiptl/real_space/logs/'
models = get_hiptl_models()

# getting command-line input
START = int(sys.argv[1])
END = int(sys.argv[2])
SNAPSHOT = int(sys.argv[3])
BOX = int(sys.argv[4])

# opening files to write to
w = hp.File(BASE+'hiptl%d_%03d.%d.%d.hdf5'%(BOX, SNAPSHOT, START, END),'w')
logfile = open(BASE+'combine%d_%03d.log'%(BOX, SNAPSHOT),'a')

# getting the filenames that we will be combining
filenos = np.arange(int(START), int(END))
files = ['hiptl%d_%03d.%d.hdf5'%(BOX, SNAPSHOT, i) for i in filenos]

logfile.write('first file: ' + files[0]+'\n')
logfile.write('last file: ' + files[-1]+'\n')

# sum each model's grid individually
for m in models:
    total = np.zeros((2048, 2048, 2048), dtype=np.float32)
    logfile.write("starting model "+m+'\n')
    logfile.write("current total sum %.4f"%(np.sum(total)))
    for i in files:
        # it is expected that the last job will have nonexistant files
        try:
            f = hp.File(BASE+i,'r')
        except IOError:
            logfile.write('did not find the file %s\n'%i)
        else:
            total += f[m][:]
            logfile.write("new sum:" + str(np.sum(total))+"\n")
            f.close()
    w.create_dataset(m, data=total, compression="gzip", compression_opts=9)
w.close()
logfile.close()
