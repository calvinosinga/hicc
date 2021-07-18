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
SAVE = '/lustre/cosinga/HI-color/results/plots/aux_plots/'

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

# removing the unresolved subhalos in both gas AND stmass
# TODO: maybe find the % of each bin that is resolved/unresolved
both_idx = lhicc.is_resolved_gas_stmass(stmass, gasmass)
stmass_both = stmass[both_idx]
gr_both = gr[both_idx]

# removing unresolved subhalos in just stellar mass
stres_idx = lhicc.is_resolved_stmass(stmass)
stmass_stres = stmass[stres_idx]
gr_stres = gr[stres_idx]

# removing unresolved subhalos in stellar mass or gas mass
res_idx = lhicc.is_resolved_nelson(stmass, gasmass)
stmass_res = stmass[res_idx]
gr_res = gr[res_idx]

# the nelson lines definitions in gr-stmass plane
fx = lambda x: 0.65 + 0.02*(x-10.28)
fxdown = lambda x: 0.625 + 0.02*(x-10.28)
fxup = lambda x: 0.675 + 0.02*(x-10.28)

# now making first histogram
stmass_both = np.log10(stmass_both)
plt.hist2d(stmass_both,gr_both,bins=50,norm=mpl.colors.LogNorm())
cbar = plt.colorbar()

# adding lines to show the definitions for blue/red
x = np.linspace(np.min(stmass_both), np.max(stmass_both))
plt.plot(x, fx(x), label='fiducial', color='red', linestyle='-')
plt.plot(x, fxdown(x), label='bluer cut', color='tomato', linestyle=':')
plt.plot(x, fxup(x), label='redder cut', color='darkred', linestyle=':')
plt.legend()
plt.xlabel('Stellar Mass')
plt.ylabel('g-r (magnitude)')
plt.title('Resolved in Both Gas and Stellar Mass')
plt.savefig(SAVE+'gas+stmass_resolved_%d_%03d.png'%(BOX,SNAPSHOT))
plt.clf()

# making second histogram
stmass_stres = np.log10(stmass_stres)
plt.hist2d(stmass_stres,gr_stres,bins=50,norm=mpl.colors.LogNorm())
cbar = plt.colorbar()

# adding lines to show the definitions for blue/red
x = np.linspace(np.min(stmass_stres), np.max(stmass_stres))
plt.plot(x, fx(x), label='0.65+0.02(M-10.28)', color='red', linestyle='-')
plt.plot(x, fxdown(x), label='0.6+0.02(M-10.28)', color='tomato', linestyle='--')
plt.plot(x, fxup(x), label='0.7+0.02(M-10.28)', color='darkred', linestyle='--')
plt.legend()
plt.xlabel('Stellar Mass')
plt.ylabel('g-r (magnitude)')
plt.title('Color Bimodality z=0 TNG100')
plt.savefig(SAVE+'stmass_resolved_%d_%03d.png'%(BOX,SNAPSHOT))
plt.clf()

# making third histogram
stmass_res = np.log10(stmass_res)
stmass_res[stmass_res==-np.inf] = 0
plt.hist2d(stmass_res,gr_res,bins=50,norm=mpl.colors.LogNorm())
cbar = plt.colorbar()

# adding lines and plot labels
x = np.linspace(np.min(stmass_res), np.max(stmass_res))
plt.plot(x, fx(x), label='0.65+0.02(M-10.28)', color='red', linestyle='-')
plt.plot(x, fxdown(x), label='0.6+0.02(M-10.28)', color='tomato', linestyle='--')
plt.plot(x, fxup(x), label='0.7+0.02(M-10.28)', color='darkred', linestyle='--')
plt.legend(loc='lower right')
plt.xlabel('Stellar Mass')
plt.ylabel('g-r (magnitude)')
plt.title('Resolved in Stellar Mass or Gas Mass')
plt.savefig(SAVE+'gas_or_stmass_resolved_%d_%03d.png'%(BOX,SNAPSHOT))
plt.clf()
