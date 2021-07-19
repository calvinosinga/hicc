from h5py._hl.files import File
from pk.pk import BOXSIZE, SNAPSHOT
import sys
import numpy as np
import h5py as hp
import matplotlib.pyplot as plt
import library_hicc.plot as lpt
from library_hicc.printer import Printer
import os
import illustris_python as il

# getting command-line arguments

FILE = sys.argv[1]
SNAPSHOT = int(sys.argv[2])
BOX = int(sys.argv[3])

# defining paths
PLOTS = '/lustre/cosinga/HI-color/results/plots/'
FIELDS = '/lustre/cosinga/HI-color/results/fields/snap_%03d/'%SNAPSHOT
LOGS = '/lustre/cosinga/HI-color/hicc/logs/slices/'
TNG = '/lustre/cosinga/tng%s'%BOX

# making printer
pnt = Printer(LOGS+"%s.log"%FILE)
pnt.write("printer made")

# getting sim data
pnt.write("getting simulation data...")
head = il.groupcat.loadHeader(TNG, SNAPSHOT)
SCALE = head['Time'] # scale factor
BOXSIZE = head['BoxSize']/1e3 * SCALE #Mpc/h - convert from comoving using a
del head

# making directory to store plots in...
if not os.path.isdir(PLOTS+"/slices/%s/"%FILE):
    os.mkdir(PLOTS+"/slices/%s/"%FILE)
    pnt.write("made directory at: "+PLOTS+"/slices/%s/"%FILE)

# getting input data
filename = FILE + "%d_%03d.final."%(BOX, SNAPSHOT)
pnt.write("opening input file at: "+FIELDS+filename+'.hdf5')
f = hp.File(FIELDS+filename+'.hdf5','r')
keylist = list(f.keys())


# now making slice plots
for key in keylist:
    field = f[key][:]
    pnt.write("making slice plot for %s"%key)
    lpt.plotslc(field, BOXSIZE, PLOTS+'/slices/%s/%s%d_%03d'%(FILE, key, BOX, SNAPSHOT))
    pnt.writeTab("successfully made the slice plot")
    count = np.count_nonzero(field)
    elms = field.shape[0]**3
    pnt.writeTab("out of %.3e elements, %.3e are occupied"%(elms, count))

f.close()
