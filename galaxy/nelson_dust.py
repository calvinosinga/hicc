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
SAVE = '/lustre/cosinga/HI-color/results/fields/snap_%03d/'%SNAPSHOT
LOG = '/lustre/cosinga/HI-color/hicc/logs/galaxy/'

# getting simulation defined constants
head = il.groupcat.loadHeader(HOME,SNAPSHOT)
LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
SCALE = head['Time'] # scale factor
BOXSIZE = head['BoxSize']/1e3 * SCALE #Mpc/h
REDSHIFT = head['Redshift']

# output files
wrs = hp.File('%snelsonrs_%s_dust%d_%03d.%dres.final.hdf5'%(SAVE,RUN,BOX,SNAPSHOT,RES), 'w')
w = hp.File('%snelson_%s_dust%d_%03d.%dres.final.hdf5'%(SAVE,RUN,BOX,SNAPSHOT, RES), 'w')
pnt = Printer(LOG+'nelson_%s_dust%d_%03d.%dres.log'%(RUN,BOX,SNAPSHOT,RES))

# input data
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

# getting dust-adjusted photometric data
dustfile = "Subhalo_StellarPhot_p07c_cf00dust_res_conv_ns1_rad30pkpc_%03d.hdf5"%SNAPSHOT
f = hp.File(HOME+"postprocessing/stellar_light/"+dustfile)
dustPhoto = f['Subhalo_StellarPhot_p07c_cf00dust_res_conv_ns1_rad30pkpc'][:]

# get the angles for each of the line-of-sights that the dust was computed for
proj = dict(dustPhoto.attrs)['projVecs']

# get the line-of-sight for the power spectrum calculation
los = np.zeros_like(proj)
los[:,AXIS] += 1

# find the projection that is the closest to the pk's line of sight
dist = np.sum((proj-los)**2,axis=1)
minidx = np.argmin(dist)

pnt.write("the projection found was "+str(proj[minidx]))

gr = dustPhoto[:, 1, minidx]-dustPhoto[:, 2, minidx]

pnt.write("writing to the redshift file %snelsonrs_%s%d_%03d.final.hdf5"%(SAVE,RUN,BOX,SNAPSHOT))
pnt.write("writing to the real space file %snelson_%s%d_%03d.final.hdf5"%(SAVE,RUN,BOX,SNAPSHOT))

counts = []
counts_names = []
def create_field(fieldname, mask):
    """
    Creates a field using the only the indices provided, saved to 
    the output file w using the fieldname as the key.
    """

    pnt.write("creating fields for %s"%fieldname)
    pnt.write("\t now creating real-space field")
    field = np.zeros(GRID, dtype=np.float32)
    
    CICW(pos[mask], field, BOXSIZE, total_mass[mask])
    pnt.write("\t now assigning the masses to the grid")
    pnt.write("\t average mass of each subhalo: %d"%np.mean(total_mass[mask]))
    w.create_dataset(fieldname, data=field, compression="gzip", compression_opts=9)
    field = np.zeros(GRID, dtype=np.float32)
    pnt.write("\t now creating redshift-space field")
    pnt.write("\t now assigning the masses to the grid")
    CICW(rspos[mask], field, BOXSIZE, total_mass[mask]) #rspos is adjusted outside of this method
    wrs.create_dataset(fieldname, data=field)
    counts.append(np.sum(mask))
    counts_names.append(fieldname)
    pnt.write("\t there are %d subhalos in this group"%np.sum(mask))
    return

pnt.write("starting subhalo grid generation")

resolved_mask = np.invert(np.isnan(gr))
red_mask = lhicc.is_red_nelson(gr,mass[:,4],RUN)
blue_mask = np.invert(red_mask)
red_mask *= resolved_mask # removes unresolved true values
blue_mask *= resolved_mask

# creating the fields
create_field("total", np.ones(total_mass.shape, dtype=bool))
create_field("resolved", resolved_mask)
create_field("unresolved", np.invert(resolved_mask))
create_field("blue", blue_mask)
create_field("red", red_mask)

pnt.write("now saving the counts")
# saving the counts
cfile = open(SAVE+"subhalo_counts%s%d_%03d.txt"%(RUN,BOX,SNAPSHOT),'w')
for c in range(len(counts)):
    cfile.write("%s %d\n"%(counts_names[c], counts[c]))

pnt.write("now closing the files")
cfile.close()
w.close()
