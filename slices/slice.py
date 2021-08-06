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
RES = int(sys.argv[4])

# defining paths
PLOTS = '/lustre/cosinga/HI-color/results/plots/'
FIELDS = '/lustre/cosinga/HI-color/results/fields/snap_%03d/'%SNAPSHOT
LOGS = '/lustre/cosinga/HI-color/hicc/logs/slices/'
TNG = '/lustre/cosinga/tng%s'%BOX

# making printer
pnt = Printer(LOGS+"%s.%dres.log"%(FILE,RES))
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
filename = FILE + "%d_%03d.%dres.final"%(BOX, SNAPSHOT,RES)
pnt.write("opening input file at: "+FIELDS+filename+'.hdf5')
f = hp.File(FIELDS+filename+'.hdf5','r')
keylist = list(f.keys())

def coarse_grid(grid):
    grid = grid[1::2,:,:]+grid[::2,:,:]
    grid = grid[:,1::2,:]+grid[:,::2,:]
    grid = grid[:,:,1::2]+grid[:,:,::2]
    return grid

# now making slice plots
for key in keylist:
    field = f[key][:]
    if not len(field.shape) == 3:
        pnt.write("field %s is not the right shape; skipping slice plot..."%key)
    else:
        count = np.count_nonzero(field)
        elms = field.shape[0]**3
        occupation = count/elms
        if occupation < .001:
            pnt.writeTab("found that this grid is too sparse, converting to a smaller grid")
            field = coarse_grid(field)
        pnt.writeTab("out of %.3e elements, %.3e are occupied"%(elms, count))
        pnt.writeTab("making slice plot for %s"%key)
        lpt.plotslc(field, BOXSIZE, PLOTS+'/slices/%s/%s%d_%03d.%dres'%(FILE, key, BOX, SNAPSHOT, RES))
        pnt.writeTab("successfully made the slice plot")


f.close()
