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
f= hp.File(HOME+'snapdir_%03d/snap_%03d.0.hdf5'%(SNAPSHOT,SNAPSHOT), 'r')
head = dict(f['Header'].attrs)

LITTLE_H = head['HubbleParam'] # 100 km/s/Mpc
SCALE = head['Time'] # scale factor
BOXSIZE = head['BoxSize']/1e3 * SCALE #Mpc/h
REDSHIFT = head['Redshift']
DMPTL = head['MassTable'][1]*1e10/LITTLE_H
# output files
wrs = hp.File('%snelsonrs_%s_ptl%d_%03d.%dres.final.hdf5'%(SAVE,RUN,BOX,SNAPSHOT,RES), 'w')
w = hp.File('%snelson_%s_ptl%d_%03d.%dres.final.hdf5'%(SAVE,RUN,BOX,SNAPSHOT,RES), 'w')
pnt = Printer(LOG+'nelson_%s_ptl%d_%03d.%dres.log'%(RUN,BOX,SNAPSHOT,RES))

# input data
flds = ['SubhaloPos', 'SubhaloStellarPhotometrics', 'SubhaloMassType', 'SubhaloVel']
sub = il.groupcat.loadSubhalos(HOME,SNAPSHOT, fields=flds)
mass = sub[flds[2]][:]*1e10/LITTLE_H # solar masses
total_mass = np.sum(mass, axis=1)
gr = sub[flds[1]][:,4] - sub[flds[1]][:,5]
pos = sub[flds[0]][:]/1e3 * SCALE # Mpc/h, 52 MB
vel = sub[flds[3]][:] # km/s, 52 MB

# pnt.write("now shifting the positions to redshift space...")
# rspos = pos_redshift_space(pos, vel, BOXSIZE, 100*LITTLE_H, REDSHIFT, AXIS)
del sub, flds

# pnt.write("writing to the redshift file %snelsonrs_%s%d_%03d.final.hdf5"%(SAVE,RUN,BOX,SNAPSHOT))
# pnt.write("writing to the real space file %snelson_%s%d_%03d.final.hdf5"%(SAVE,RUN,BOX,SNAPSHOT))
counts = []
counts_names = []
# def create_field(fieldname, mask):
#     """
#     Creates a field using the only the indices provided, saved to 
#     the output file w using the fieldname as the key.
#     """
#     pnt.write("creating fields for %s"%fieldname)
#     pnt.write("\t now creating real-space field")
#     field = np.zeros(GRID, dtype=np.float32)
    
#     CICW(pos[mask], field, BOXSIZE, total_mass[mask])

#     pnt.write("\t now assigning the masses to the grid")
#     pnt.write("\t average mass of each subhalo: %d"%np.mean(total_mass[mask]))

#     w.create_dataset(fieldname, data=field, compression="gzip", compression_opts=9)
#     field = np.zeros(GRID, dtype=np.float32)
#     pnt.write("\t now creating redshift-space field")
#     pnt.write("\t now assigning the masses to the grid")
#     CICW(rspos[mask], field, BOXSIZE, total_mass[mask]) #rspos is adjusted outside of this method
#     wrs.create_dataset(fieldname, data=field)
#     counts.append(np.sum(mask))
#     counts_names.append(fieldname)
#     pnt.write("\t there are %d subhalos in this group"%np.sum(mask))
#     return

pnt.write("starting subhalo grid generation")

resolved_mask = lhicc.is_resolved_stmass(mass[:,4])
red_mask = lhicc.is_red_nelson(gr,mass[:,4],RUN)
blue_mask = np.invert(red_mask)
red_mask *= resolved_mask # removes unresolved true values
blue_mask *= resolved_mask

# creating the fields
# create_field("total", np.ones(total_mass.shape, dtype=bool))
# create_field("resolved", resolved_mask)
# create_field("unresolved", np.invert(resolved_mask))
# create_field("blue", blue_mask)
# create_field("red", red_mask)

# creating the fields for blue-ptl, red-ptl, resolved-ptl
masks = [blue_mask, red_mask, resolved_mask]
names = ['blue', 'red', 'resolved']
ptltypes = [0,1,4,5]
for m in range(len(masks)):
    grid = np.zeros(GRID, dtype=np.float32)
    gridrs = np.zeros(GRID, dtype=np.float32)
    for i in range(len(masks[m])):
        if masks[m][i]:
            for p in ptltypes:
                print(p)
                
                if p==1:
                    ptldata=il.snapshot.loadSubhalo(HOME, SNAPSHOT, i, p, fields=['Coordinates','Velocities'])
                    mass = np.ones(len(ptldata['Coordinates']), dtype=np.float32)*DMPTL
                else:
                    ptldata = il.snapshot.loadSubhalo(HOME, SNAPSHOT, i, p, fields=['Masses','Coordinates','Velocities'])
                    mass = ptldata['Masses']*1e10/LITTLE_H
                pos = ptldata['Coordinates']/1e3 * SCALE
                vel = ptldata['Velocities'] * np.sqrt(SCALE)
                
                CICW(pos, grid, BOXSIZE, mass)

                rspos = pos_redshift_space(pos, vel, BOXSIZE, LITTLE_H*100, REDSHIFT, AXIS)
                CICW(rspos, gridrs, BOXSIZE, mass)

    w.create_dataset(names[m], data=grid, compression="gzip", compression_opts=9)
    wrs.create_dataset(names[m], data=gridrs, compression="gzip", compression_opts=9)



# pnt.write("now saving the counts")
# # saving the counts
# cfile = open(SAVE+"subhalo_counts%s%d_%03d.txt"%(RUN,BOX,SNAPSHOT),'w')
# for c in range(len(counts)):
#     cfile.write("%s %d\n"%(counts_names[c], counts[c]))

# pnt.write("now closing the files")
# cfile.close()
w.close()
wrs.close()
