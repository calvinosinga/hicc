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
from library_hicc.redshift_space import pos_redshift_space
import library_hicc.color as lhicc
import illustris_python as il
from library_hicc.printer import Printer

# reading command line inputs
SNAPSHOT = int(sys.argv[1])
BOX = int(sys.argv[2])
RUN = sys.argv[3] # low,mid,high to test sensitivity to color definition
RES = int(sys.argv[4]) # resolution of the grid
AXIS = int(sys.argv[5])
GRID = (RES,RES,RES)

# defining needed paths
HOME = '/lustre/cosinga/tng%d/'%BOX
SAVE = '/lustre/cosinga/final_fields/'
LOG = '/lustre/cosinga/hicc/logs/'

# getting simulation defined constants
head = il.groupcat.loadHeader(HOME,SNAPSHOT)
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
SCALE = head['Time'] # scale factor
BOXSIZE = head['BoxSize']/1e3 * SCALE #Mpc/h
REDSHIFT = head['Redshift']

# output files
wrs = hp.File('%snelsonrs_%s%d_%03d.final.hdf5'%(SAVE,RUN,BOX,SNAPSHOT), 'w')
w = hp.File('%snelson_%s%d_%03d.final.hdf5'%(SAVE,RUN,BOX,SNAPSHOT), 'w')
pnt = Printer(LOG+'nelson_%s%d_%03d.log'%(RUN,BOX,SNAPSHOT))
# input data
pnt.write("starting subhalo grid generation")
flds = ['SubhaloPos', 'SubhaloStellarPhotometrics', 'SubhaloMassType', 'SubhaloVel']
sub = il.groupcat.loadSubhalos(HOME,SNAPSHOT, fields=flds)
mass = sub[flds[2]][:]*1e10/LITTLE_H # solar masses
total_mass = np.sum(mass, axis=1)
gr = sub[flds[1]][:,4] - sub[flds[1]][:,5]
pos = sub[flds[0]][:]/1e3 * SCALE # Mpc/h, 52 MB
vel = sub[flds[3]][:] # km/s, 52 MB
pnt.write("now shifting the positions to redshift space...")
rspos = pos_redshift_space(pos, vel, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)
del sub, flds



pnt.write("writing to the redshift file %snelsonrs_%s%d_%03d.final.hdf5"%(SAVE,RUN,BOX,SNAPSHOT))
pnt.write("writing to the real space file %snelson_%s%d_%03d.final.hdf5"%(SAVE,RUN,BOX,SNAPSHOT))
counts = []
counts_names = []
def create_field(fieldname, idx):
    """
    Creates a field using the only the indices provided, saved to 
    the output file w using the fieldname as the key.
    """
    def print_mem(name, obj):
        mem = sys.getsizeof(obj)
        mem = mem/1e6
        pnt.write("\t %s is %.3eMB"%(name, mem))
        return
    pnt.write("creating fields for %s"%fieldname)
    pnt.write("\t now creating real-space field")
    field = np.zeros(GRID, dtype=np.float32)
    
    CICW(pos[idx], field, BOXSIZE, total_mass[idx])
    pnt.write("\t now assigning the masses to the grid")
    pnt.write("\t average mass of each subhalo: %d"%np.mean(total_mass[idx]))
    print_mem("mass", total_mass[idx])
    print_mem("field", field)
    print_mem("position", pos[idx])
    w.create_dataset(fieldname, data=field, compression="gzip", compression_opts=9)
    field = np.zeros(GRID, dtype=np.float32)
    pnt.write("\t now creating redshift-space field")
    pnt.write("\t now assigning the masses to the grid")
    CICW(rspos[idx], field, BOXSIZE, total_mass[idx]) #rspos is adjusted outside of this method
    wrs.create_dataset(fieldname, data=field)
    counts.append(np.sum(idx))
    counts_names.append(fieldname)
    pnt.write("\t there are %d subhalos in this group"%np.sum(idx))
    return

resolved_idx = lhicc.is_resolved_stmass(mass[:,4])
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

pnt.write("now saving the counts")
# saving the counts
cfile = open(SAVE+"subhalo_counts%s%d_%03d.txt"%(RUN,BOX,SNAPSHOT),'w')
for c in range(len(counts)):
    cfile.write("%s %d\n"%(counts_names[c], counts[c]))

pnt.write("now closing the files")
cfile.close()
w.close()
