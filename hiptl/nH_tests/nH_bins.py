from hisubhalo.hisubhalo import SCALE
import numpy as np
import h5py as hp
import sys
from library_hicc.mas import CICW
from library_hicc.redshift_space import pos_redshift_space
from library_hicc.models import get_hiptl_models

# reading command line inputs
CHUNK = int(sys.argv[1])
SNAPSHOT = int(sys.argv[2])
BOX = int(sys.argv[3])
RES = int(sys.argv[4]) # resolution of the grid
AXIS = int(sys.argv[5]) # if -1, not in redshift space
IN_RS_SPACE = not (AXIS == -1)

# defining needed paths
PTLPATH = '/lustre/cosinga/tng%d/snapdir_%03d/'%(BOX, SNAPSHOT) # where the ptl files are saved
OUTPATH = '/lustre/cosinga/hiptl_output/' # where to save the output
HIPATH = '/lustre/diemer/illustris/hih2/'  # where the hiptl files are saved

# input files
ptlfile = hp.File(PTLPATH+"snap_%03d.%d.hdf5" %(SNAPSHOT, CHUNK), 'r')
hih2file = hp.File(HIPATH+"hih2_particles_%03d.%d.hdf5" %(SNAPSHOT, CHUNK), 'r')

# output files
if IN_RS_SPACE:
    w = hp.File(OUTPATH+'nHrs%d_%03d.%d.hdf5' %(BOX, SNAPSHOT, CHUNK), 'w')
else:
    w = hp.File(OUTPATH+'nH%d_%03d.%d.hdf5' %(BOX, SNAPSHOT, CHUNK), 'w')

# getting author-defined constants (these COULD change but are not expected to)
GRID = (RES,RES,RES)

# getting the needed simulation-defined constants
head = dict(ptlfile['Header'].attrs)
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
BOXSIZE = head['BoxSize']/1e3 #Mpc/h
REDSHIFT = head['Redshift']
SCALE_FACTOR = head['Time']
DMPTL = head['MassTable'][1] *1e10/LITTLE_H # solar masses
ptltype = [0,1,4,5]
ptltype = ['PartType'+str(i) for i in ptltype]
nptl = head['NumPart_ThisFile']

# testing the different hydrogen densities

dendec =[0,1e-4,1e-3,1e-2,1e-1,1,10,100,np.inf]
counts = np.zeros_like(dendec)
p = 'PartType0'
mass = ptlfile[p]['Masses'][:]*1e10/LITTLE_H # solar masses
density = ptlfile[p]['Density'][:]*1e10/LITTLE_H/SCALE_FACTOR**3 # solar masses per kpc/h cubed
# tothyd = ptlfile[p]['GFM_Metals'][:,0] # fraction of cell that is hydrogen
f_neut_h = hih2file['PartType0']['f_neutral_H'][:]
# replacing the cells with negative fractions with zero
f_neut_h = np.where(f_neut_h>=0, f_neut_h, np.zeros_like(f_neut_h))
kpctocm = 3.086e21
smtog = 1.989e33
m_p=1.673e-24
fac = 1/kpctocm**3*LITTLE_H**3*smtog/m_p

n_h = density*f_neut_h*fac
pos = ptlfile[p]['Coordinates'][:]/1e3 # Mpc/h
vel = ptlfile[p]['Velocities'][:] * np.sqrt(SCALE_FACTOR) # km/s
models = get_hiptl_models()
models.append("total")
print("now binning according to number density, there are " + str(n_h.shape)+ " cells.")
for m in models:
    if not m is "total":
        molfrac = hih2file[p][m][:]

    for d in range(len(dendec)-1):
        field = np.zeros(GRID, dtype=np.float32)
        lo = dendec[d]
        hi = dendec[d+1]
        print("getting low mask: %.3e"%lo)
        mask1 = n_h >= lo
        print("low mask has %d cells in it"%np.sum(mask1))
        print("getting high mask: %.3e"%hi)
        mask2 = n_h < hi
        print("high mask has %d cells in it"%np.sum(mask2))
        mask = mask1 & mask2
        counts[d] = np.sum(mask)
        print("in the <%.3e bin there are %d cells, has average %.4e"%(hi, counts[d], np.mean(n_h[mask])))

        posmask = pos[mask,:]
        if m is "total":
            massmask = mass[mask] * f_neut_h[mask]
        else:
            massmask = mass[mask] * f_neut_h[mask] * (1-molfrac[mask])
        if IN_RS_SPACE:
            velmask = vel[mask,:]
            posmask = pos_redshift_space(posmask, velmask, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)
        CICW(posmask,field,BOXSIZE,massmask)
        w.create_dataset("%s_<%.3e"%(m,dendec[d]), data=field, compression="gzip", compression_opts=9)
w.create_dataset("bin_counts",data=counts)
w.close()
