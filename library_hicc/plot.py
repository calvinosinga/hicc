#!/usr/bin/env python3


import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import copy
import h5py as hp

"""
This is meant to create and store plots of the pks and slices of their fields to check that they make sense.
"""

def plot1Dpk(k, pk, res, boxsize, plotname):
    
    # the nyquist frequency is N*pi/L
    nyq = res*np.pi/boxsize
    maxy = np.max(pk)
    plt.plot(k,pk)
    plt.loglog()
    plt.ylabel('P(k) (h/Mpc)^3')
    plt.xlabel('k (h/Mpc)')
    plt.plot((nyq,nyq), (0,maxy), 'k:')
    plt.savefig(plotname+'.png')
    plt.clf()
    return

def plot2Dpk(kpar, kper, pk, plotname):
    kpar = np.unique(kpar)
    kper = np.unique(kper)
    levs = [-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4]
    pk = np.log10(np.reshape(pk,(len(kper),len(kpar))))
    cmap = copy.copy(mpl.cm.get_cmap("plasma"))
    cmap.set_under('w')
    paridx = np.where(kpar>4)[0][0]
    peridx = np.where(kper>4)[0][0]
    KPAR, KPER = np.meshgrid(kpar[:paridx], kper[:peridx])
    plt.contour(KPAR,KPER, pk[:paridx, :peridx], vmin=-2, vmax=4, levels=levs, colors='k', linestyles='solid')
    plt.imshow(pk[:paridx,:peridx], vmin=-2, vmax=4, extent=(0,kpar[paridx-1],0,kper[peridx-1]), origin='lower', cmap=cmap)
    plt.xlabel("kpar (h/Mpc)")
    plt.ylabel("kper (h/Mpc)")
    plt.colorbar()
    plt.savefig(plotname+'.png')
    plt.clf()
    return

def plotslc(grid, boxsize, plotname):
    cmap = copy.copy(mpl.cm.get_cmap("plasma"))
    # cmap.set_under('w')
    dim = grid.shape[0]
    slcidx = int(0.1*dim)
    mid = int(dim/2)
    slc = np.log10(np.sum(grid[:, mid-slcidx:mid+slcidx, :], axis=1))
    w = hp.File(plotname+'.hdf5', 'w')
    w.create_dataset('slice', data=slc)
    w.close()
    plt.imshow(slc, extent=(0,boxsize, 0, boxsize), origin='lower', cmap=cmap)
    plt.xlabel("x (Mpc/h)")
    plt.ylabel("y (Mpc/h)")
    cbar = plt.colorbar()
    cbar.set_label("Mass (Solar Masses)")
    plt.savefig(plotname+'.png')
    plt.clf()
    return
    
def plotCorr(r, xi, res, boxsize, plotname):
    nyq = 2*boxsize/res
    maxy = np.max(xi)
    plt.plot(r, xi)
    plt.xlabel("r (Mpc/h)")
    plt.ylabel("xi (Mpc/h)")
    plt.loglog()
    plt.plot((nyq, nyq), (0,maxy), 'k:')
    plt.savefig(plotname+'.png')
    plt.clf()
    return
