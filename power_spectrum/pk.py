#!/usr/bin/env python3
"""
Calculates the power spectrum of the fields given by a text file.
The name of the text file is given as command-line input, .
"""

# import statements
import sys
import numpy as np
import h5py as hp
import illustris_python as il
from Pk_library import Pk
from Pk_library import Xpk

# reading command line inputs
PKPATH = sys.argv[1]
LINENO = int(sys.argv[2])
SAVENAME = sys.argv[3]
SNAPSHOT = int(sys.argv[4])
BOX = sys.argv[5]

# defining basic paths
HOME = '/lustre/cosinga/final_fields/'
TNG = '/lustre/cosinga/tng%d'%BOX

# getting author-defined constants
MAS = 'CIC'

# getting simulation defined constants
head = il.groupcat.loadHeader(TNG, SNAPSHOT)
BOXSIZE = head['BoxSize']/1e3 # Mpc/h
del head

# getting input files
pklist = list(open(HOME+PKPATH, 'r'))
pklist = pklist[LINENO][:-1].split()

# interpretting input
if pklist[0] == 'cross':
    IS_XPK = True
    IN_RS = False
    pk1 = pklist[0]
    key1 = pklist[1]
    pk2 = pklist[2]
    key2 = pklist[3]
    dim = int(pklist[4])
    axis = int(pklist[5])
elif pklist[0] == 'auto':
    IS_XPK = False
    IN_RS = False
    pk1 = pklist[0]
    key1 = pklist[1]
    dim = int(pklist[2])
    axis = int(pklist[3])
else:
    raise ValueError("incorrect first input on line %d"%LINENO)

del pklist

def to_overdensity(field):
    field = field/np.mean(field).astype(np.float32)
    field=field - 1

if IS_XPK:
    f1 = hp.File(HOME+pk1,'r')
    f2 = hp.File(HOME+pk2,'r')
    field1 = f1[key1][:]
    field2 = f2[key2][:]
    to_overdensity(field1)
    to_overdensity(field2)
    res = XPk([field1,field2], BOXSIZE, axis=axis, MAS=[MAS, MAS])
    if dim == 1:
        xpk = np.transpose([res.k3D, res.XPk[:,0,0]])
        np.savetxt('%s/pk/%s%d.txt'%(HOME, SAVENAME, LINENO), xpk)
    elif dim == 2:
        x2pk = np.transpose([res.kpar, res.kper, res.PkX2D[:,0,0]])
        np.savetxt('%s/2dpk/%s%d.txt'%(HOME, SAVENAME, LINENO), x2pk)
    else:
        xpk = np.transpose([res.k3D, res.XPk[:,0,0]])
        np.savetxt('%s/pk/%s%d.txt'%(HOME, SAVENAME, LINENO), xpk)
        x2pk = np.transpose([res.kpar, res.kper, res.PkX2D[:,0,0]])
        np.savetxt('%s/2dpk/%s%d.txt'%(HOME, SAVENAME, LINENO), x2pk)
    f1.close()
    f2.close()
else:
    f1 = hp.File(HOME+pk1,'r')
    field1 = f1[key1][:]
    to_overdensity(field1)
    res = Pk(field1, BOXSIZE, axis=axis, MAS=MAS)
    if dim == 1:
        pk = np.transpose([res.k3D, res.Pk[:,0]])
        np.savetxt('%s/pk/%s%d.txt'%(HOME, SAVENAME, LINENO), xpk)
    elif dim == 2:
        d2pk = np.transpose([res.kpar, res.kper, res.Pk2D[:,0]])
        np.savetxt('%s/2dpk/%s%d.txt'%(HOME, SAVENAME, LINENO), x2pk)
    else:
        xpk = np.transpose([res.k3D, res.Pk[:,0]])
        np.savetxt('%s/pk/%s%d.txt'%(HOME, SAVENAME, LINENO), xpk)
        x2pk = np.transpose([res.kpar, res.kper, res.Pk2D[:,0]])
        np.savetxt('%s/2dpk/%s%d.txt'%(HOME, SAVENAME, LINENO), x2pk)
    f1.close()
