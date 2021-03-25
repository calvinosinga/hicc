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
from Pk_library import XPk

# reading command line inputs
AUTO_OR_XPK = sys.argv[1]
SNAPSHOT = int(sys.argv[2])
BOX = int(sys.argv[3])
DIM = int(sys.argv[4])
AXIS = int(sys.argv[5])
FILE1 = sys.argv[6]
if AUTO_OR_XPK == "cross":
    FILE2 = sys.argv[7]
    IS_XPK = True
elif AUTO_OR_XPK == "auto":
    IS_XPK = False
else:
    raise ValueError("incorrect input: first arg should be auto or cross")


# defining basic paths
HOME = '/lustre/cosinga/final_fields/'
TNG = '/lustre/cosinga/tng%d'%BOX

# getting author-defined constants
MAS = 'CIC'

# getting simulation defined constants
head = il.groupcat.loadHeader(TNG, SNAPSHOT)
BOXSIZE = head['BoxSize']/1e3 # Mpc/h
del head
print("The boxsize is %d"%BOXSIZE)

def to_overdensity(field):
    field = field/BOXSIZE**3
    field = field/np.mean(field).astype(np.float32)
    field=field - 1
    if not field.dtype == np.float32:
        print("turning into float32s")
        field.astype(np.float32)
    

if IS_XPK:
    print("calculating the cross-power for %s, %s"%(FILE1, FILE2))
    # output file
    w = hp.File(HOME+'pk/%s-%s%d_%03d.%dDxpk.hdf5'%(FILE1,FILE2,BOX,SNAPSHOT,DIM),'w')

    # input files, loop over all of the keys in each
    f1 = hp.File(HOME+FILE1+'%d_%03d.final.hdf5'%(BOX, SNAPSHOT),'r')
    key1list = list(f1.keys())
    f2 = hp.File(HOME+FILE2+'%d_%03d.final.hdf5'%(BOX, SNAPSHOT),'r')
    key2list = list(f2.keys())
    for key1 in key1list:
        for key2 in key2list:
            print("starting xpk calculation for %s, %s"%(key1, key2))
            # get the fields
            field1 = f1[key1][:]
            field2 = f2[key2][:]
            
            # convert them to overdensities
            to_overdensity(field1)
            to_overdensity(field2)
            # compute the xpk
            res = XPk([field1,field2], BOXSIZE, axis=AXIS, MAS=[MAS, MAS])

            # if this is the first calculation, save the wavenumbers
            if key1 == key1list[0] and key2 == key2list[0]:
                # if in 1D, just need 1 k, otherwise need kper and kpar
                if DIM == 1:
                    w.create_dataset("k",data=res.k3D)
                elif DIM == 2:
                    w.create_dataset("kper", data=res.kper)
                    w.create_dataset("kpar", data=res.kpar)
            
            # save xpk result    
            if DIM == 1:
                w.create_dataset("%s-%s"%(key1, key2),data=res.XPk[:,0,0])
            elif DIM == 2:
                w.create_dataset("%s-%s"%(key1, key2), data=res.PkX2D[:,0,0])
    f1.close()
    f2.close()
    w.close()
else:# auto power spectrum
    print("starting procedure for the auto power for %s"%FILE1)
    # output file
    w = hp.File(HOME+'pk/%s%d_%03d.%dDpk.hdf5'%(FILE1,BOX,SNAPSHOT,DIM),'w')

    # input data
    f1 = hp.File(HOME+FILE1+'%d_%03d.final.hdf5'%(BOX, SNAPSHOT),'r')
    keylist = list(f1.keys())
    
    # iterate over each k, save its auto power to file
    for k in keylist:
        print("calculating pk for %s"%k)
        field1 = f1[k][:]
        to_overdensity(field1)
        res = Pk(field1, BOXSIZE, axis=AXIS, MAS=MAS)

        # if first calculation, save the wavenumbers
        if k == keylist[0]:
            # if in 1D, just need 1 k, otherwise need kper and kpar
            if DIM == 1:
                w.create_dataset("k",data=res.k3D)
            elif DIM == 2:
                w.create_dataset("kper", data=res.kper)
                w.create_dataset("kpar", data=res.kpar)
        
        # save pk result
        if DIM == 1:
            w.create_dataset(k,data= res.Pk[:,0])
        elif DIM == 2:
            w.create_dataset(k,data=res.Pk2D[:])
    f1.close()
    w.close()
