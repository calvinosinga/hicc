#!/usr/bin/env python3
"""
This file makes color-stellar mass plots
"""
# import statements
import illustris_python as il
import numpy as np
import matplotlib.pyplot as plt
import sys
import library_hicc.color as lhicc
import matplotlib as mpl

# reading command-line inputs
SNAPSHOT = int(sys.argv[1])
BOX = int(sys.argv[2])

# defining needed paths
HOME = '/lustre/cosinga/tng%d/'%BOX
SAVE = '/lustre/cosinga/final_fields/plots/'

# getting input files
fields = ['SubhaloMassType','SubhaloStellarPhotometrics']
head = il.groupcat.loadHeader(HOME, SNAPSHOT)
sub = il.groupcat.loadSubhalos(HOME, SNAPSHOT, fields=fields)

# simulation-defined constants
LITTLE_H = head['HubbleParam']

# getting needed data to make plots
stmass = sub[fields[0]][:,4]*1e10/LITTLE_H
gasmass = sub[fields[0]][:,0]*1e10/LITTLE_H
gr = sub[fields[1]][:,4] - sub[fields[1]][:,5]

# removing the unresolved subhalos
# TODO: maybe find the % of each bin that is resolved/unresolved
res_idx = lhicc.is_resolved_nelson(stmass, gasmass)
stmass = stmass[res_idx]
gasmass = gasmass[res_idx]
gr = gr[res_idx]

# the nelson lines definitions in gr-stmass plane
fx = lambda x: 0.65 + 0.02*(x-10.28)
fxdown = lambda x: 0.625 + 0.02*(x-10.28)
fxup = lambda x: 0.675 + 0.02*(x-10.28)

# now making histogram
plt.hist2d(stmass,gr,bins=50,norm=mpl.colors.LogNorm())
cbar = plt.colorbar()

# adding lines to show the definitions for blue/red
x = np.linspace(np.min(stmass), np.max(stmass))
plt.plot(x, fx(x), label='fiducial', color='red', linestyle='-')
plt.plot(x, fxdown(x), label='bluer cut', color='tomato', linestyle=':')
plt.plot(x, fxup(x), label='redder cut', color='darkred', linestyle=':')
plt.legend()
plt.xlabel('Stellar Mass')
plt.ylabel('g-r (magnitude)')
plt.title('Blue and Red Subhalo Populations')
plt.savefig(SAVE+'gr-stmass_%d_%03d.png')