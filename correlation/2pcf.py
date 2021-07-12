#!/usr/bin/env python3

import sys
import numpy as np
import h5py as hp
import time
import illustris_python as il
from Pk_library import Xi
from Pk_library import XXi

# reading command-line inputs
AUTO_OR_CROSS = sys.argv[1]
SNAPSHOT = int(sys.argv[2])
BOX = int(sys.argv[3])
AXIS = int(sys.argv[4])
FILE1 = sys.argv[5]
IS_CROSS = AUTO_OR_CROSS == "cross"
if IS_CROSS:
    FILE2 = sys.argv[6]
    SAME_FILE = FILE1 == FILE2

# defining basic paths
HOME = '/lustre/cosinga/final_fields/'
TNG = '/lustre/cosinga/tng%d'%BOX

head = il.groupcat.loadHeader(TNG, SNAPSHOT)
SCALE = head['Time'] # scale factor
BOXSIZE = head['BoxSize']/1e3 * SCALE #Mpc/h
del head
MAS = "CIC"
def to_overdensity(field):
    field = field/BOXSIZE**3
    field = field/np.mean(field).astype(np.float32)
    field = field - 1
    if not field.dtype == np.float32:
        print("field not float32s... changing dtype back to float32")
        field.astype(np.float32)
    return field

if not IS_CROSS:
    # calculating auto correlation
    print("calculating auto correlation for %s"%FILE1)
    # creating output file
    w = hp.File(HOME+'pk/%s%d_%03d.corr.hdf5'%(FILE1,BOX,SNAPSHOT),'w')

    # getting input data
    f = hp.File(HOME+FILE1+'%d_%03d.final.hdf5'%(BOX, SNAPSHOT),'r')
    keylist = list(f.keys())

    for k in keylist:
        field = f[k][:]
        if not len(field.shape) == 3:
            print("skipping corr calc, field %s is not the right shape"%k)
        else:
            print("calculating correlation for %s"%k)
            field = to_overdensity(field)
            res = Xi(field, BOXSIZE, axis = AXIS, MAS=MAS)

            if k == keylist[0]:
                print("saving r...")
                w.create_dataset("r",data=res.r3D)
            w.create_dataset(k, data=res.xi[:,0])
    
    w.close()
    f.close()


else:
    # calculating cross-correlation
    print("calculating X-correlation for %s, %s"%(FILE1, FILE2))
    # creating output file
    w = hp.File(HOME+'pk/%s-%s%d_%03d.Xcorr.hdf5'%(FILE1,FILE2,BOX,SNAPSHOT),'w')

    # getting input data
    f1 = hp.File(HOME+FILE1+'%d_%03d.final.hdf5'%(BOX, SNAPSHOT),'r')
    f2 = hp.File(HOME+FILE2+'%d_%03d.final.hdf5'%(BOX, SNAPSHOT),'r')
    keylist1 = list(f1.keys())
    keylist2 = list(f2.keys())

    for k1 in keylist1:
        for k2 in keylist2:
            field1 = f1[k1][:]
            field2 = f2[k2][:]
            if not (len(field1.shape) == 3 or len(field2.shape) == 3):
                print("skipping corr calc, field %s-%s is not the right shape"%(k1,k2))
            elif SAME_FILE and k1==k2:
                print("skipping corr calc for %s: would be auto correlation..."%k1)
            else:
                print("calculating correlation for %s-%s"%(k1,k2))
                field1 = to_overdensity(field1)
                field2 = to_overdensity(field2)
                res = XXi(field1, field2, BOXSIZE, axis = AXIS, MAS=MAS)

                if k1 == keylist1[0] and k2 == keylist2[0]:
                    print("saving r...")
                    w.create_dataset("r",data=res.r3D)
                w.create_dataset("%s-%s"%(k1,k2), data=res.xi[:,0])
    
    w.close()
    f1.close()
    f2.close()



