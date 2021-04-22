#!/usr/bin/env python3
"""
Creates 2048^3 grids using the subhalo catalogue, separating them into red and blue populations.
"""

# import statements
import numpy as np
import h5py as hp
import sys
# import MAS_library as masl
from library_hicc.mas import CICW
import redshift_space_library as rsl
import library_hicc.color as lhicc
import illustris_python as il

# reading command line inputs, a fourth gives axis of redshift space
SNAPSHOT = int(sys.argv[1])
BOX = int(sys.argv[2])
RUN = sys.argv[3]
if len(sys.argv) > 4:
    AXIS = int(sys.argv[4])
    IN_RS_SPACE = True
else:
    IN_RS_SPACE = False

# defining needed paths
HOME = '/lustre/cosinga/tng%d/'%BOX
SAVE = '/lustre/cosinga/final_fields/'

# assigning author-defined constants (not expected to change)
GRID = (2048,2048,2048)
MAS = 'CIC'

# getting simulation defined constants
head = il.groupcat.loadHeader(HOME,SNAPSHOT)
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
BOXSIZE = head['BoxSize']/1e3 #Mpc/h
REDSHIFT = head['Redshift']

# input data
flds = ['SubhaloCM', 'SubhaloStellarPhotometrics', 'SubhaloMassType', 'SubhaloVel']
sub = il.groupcat.loadSubhalos(HOME,SNAPSHOT, fields=flds)
mass = sub[flds[2]][:]*1e10/LITTLE_H # solar masses
total_mass = np.sum(mass, axis=1)
gr = sub[flds[1]][:,4] - sub[flds[1]][:,5]
pos = sub[flds[0]][:]/1e3 # Mpc/h
vel = sub[flds[3]][:] # km/s
del sub, flds

# if we are working in redshift-space, shift the positions using the velocities
# then create the output file, so the names are different
if IN_RS_SPACE:
    rsl.pos_redshift_space(pos, vel, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)
    w = hp.File('%snelsonrs_%s%d_%03d.final.hdf5'%(SAVE,RUN,BOX,SNAPSHOT), 'w')
else:
    w = hp.File('%snelson_%s%d_%03d.final.hdf5'%(SAVE,RUN,BOX,SNAPSHOT), 'w')
del vel

counts = []
counts_names = []
def create_field(fieldname, idx):
    """
    Creates a field using the only the indices provided, saved to 
    the output file w using the fieldname as the key.
    """
    print("creating field for %s"%fieldname)
    field = np.zeros(GRID, dtype=np.float32)
    CICW(pos[idx], field, BOXSIZE, total_mass[idx])
    w.create_dataset(fieldname, data=field, compression="gzip", compression_opts=9)
    counts.append(np.sum(idx))
    counts_names.append(fieldname)
    return

resolved_idx = lhicc.is_resolved_stmass(mass[:,4],mass[:,0])
red_idx = lhicc.is_red_nelson(gr,mass[:,4],RUN)
blue_idx = np.invert(red_idx)
red_idx *= resolved_idx # removes unresolved true values
blue_idx *= resolved_idx

# creating the fields
create_field("total", np.ones(total_mass.shape, dtype=bool))
create_field("resolved", resolved_idx)
create_field("unresolved", np.invert(resolved_idx))
create_field("blue", blue_idx)
create_field("red", red_idx)

# saving the counts
cfile = open(SAVE+"subhalo_counts%d_%03d.txt"%(BOX,SNAPSHOT),'w')
for c in range(len(counts)):
    cfile.write("%s %d\n"%(counts_names[c], counts[c]))
cfile.close()
w.close()
