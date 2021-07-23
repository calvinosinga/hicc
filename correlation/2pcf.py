#!/usr/bin/env python3

import sys
import numpy as np
import h5py as hp
import time
import illustris_python as il
from Pk_library import Xi
from Pk_library import XXi
from library_hicc import printer
import os
import library_hicc.plot as lpt

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
SAVE = '/lustre/cosinga/HI-color/results/corr/'
HOME = '/lustre/cosinga/HI-color/results/fields/snap_%03d/'%SNAPSHOT
TNG = '/lustre/cosinga/tng%d'%BOX
LOG = '/lustre/cosinga/HI-color/hicc/logs/corr/'
PLOTS = '/lustre/cosinga/HI-color/results/plots/corr/'

# getting simulation data
head = il.groupcat.loadHeader(TNG, SNAPSHOT)
SCALE = head['Time'] # scale factor
BOXSIZE = head['BoxSize']/1e3 * SCALE #Mpc/h
del head

# author-given mass-assignment scheme
MAS = "CIC"

# creating a directory to store plots if necessary

def to_overdensity(field, pnt):
    pnt.writeTab("calculating overdensity...")
    field = field/BOXSIZE**3
    field = field/np.mean(field).astype(np.float32)
    field = field - 1
    if not field.dtype == np.float32:
        pnt.writeTab("field not float32s... changing dtype back to float32")
        field.astype(np.float32)
    return field

if not IS_CROSS:
    # creating output files
    w = hp.File(SAVE+'auto/%s%d_%03d.corr.hdf5'%(FILE1,BOX,SNAPSHOT),'w')
    pnt = printer.Printer(LOG+'%s%d_%03d.corr.log'%(FILE1,BOX,SNAPSHOT))

    # calculating auto correlation
    pnt.write("calculating auto correlation for %s"%FILE1)
    
    # getting input data
    f = hp.File(HOME+FILE1+'%d_%03d.final.hdf5'%(BOX, SNAPSHOT),'r')
    keylist = list(f.keys())

    # creating a directory to save plots to...
    if not os.path.isdir(PLOTS+'%s/'%FILE1):
        pnt.write("creating corr plot directory at: "+PLOTS+'1Dpk/%s/'%FILE1)
        os.mkdir(PLOTS+'%s/'%FILE1)
    
    
    for k in keylist:
        field = f[k][:]
        if not len(field.shape) == 3:
            pnt.write("skipping corr calc, field %s is not the right shape"%k)
        else:
            pnt.write("calculating correlation for %s"%k)
            field = to_overdensity(field, pnt)
            res = Xi(field, BOXSIZE, axis = AXIS, MAS=MAS)

            if k == keylist[0]:
                pnt.writeTab("saving r...")
                w.create_dataset("r",data=res.r3D)
            pnt.writeTab("saving xi result...")
            w.create_dataset(k, data=res.xi[:,0])

            pnt.writeTab("making correlation plot...")
            plotname = PLOTS+'%s/%s%d_%03d'%(FILE1, k, BOX, SNAPSHOT)
            lpt.plotCorr(res.r3D, res.xi[:,0], field.shape[0], BOXSIZE, plotname)
    
    w.close()
    f.close()


else:
    # creating output files
    w = hp.File(SAVE+'cross/%s-%s%d_%03d.Xcorr.hdf5'%(FILE1,FILE2,BOX,SNAPSHOT),'w')
    pnt = printer.Printer(LOG+'%s-%s%d_%03d.corr.log'%(FILE1,FILE2,BOX,SNAPSHOT))

    # calculating cross-correlation
    pnt.write("calculating X-correlation for %s, %s"%(FILE1, FILE2))

    # getting input data
    f1 = hp.File(HOME+FILE1+'%d_%03d.final.hdf5'%(BOX, SNAPSHOT),'r')
    f2 = hp.File(HOME+FILE2+'%d_%03d.final.hdf5'%(BOX, SNAPSHOT),'r')
    keylist1 = list(f1.keys())
    keylist2 = list(f2.keys())

    # making a directory to save the plots to...
    if not os.path.isdir(PLOTS+'%s-%s/'%(FILE1, FILE2)):
        os.mkdir(PLOTS+'%s-%s/'%(FILE1, FILE2))
    

    for k1 in keylist1:
        for k2 in keylist2:
            field1 = f1[k1][:]
            field2 = f2[k2][:]
            if not (len(field1.shape) == 3 or len(field2.shape) == 3):
                pnt.write("skipping corr calc, field %s-%s is not the right shape"%(k1,k2))
            elif SAME_FILE and k1==k2:
                pnt.write("skipping corr calc for %s-%s: would be auto correlation..."%(k1,k2))
            else:
                field1 = to_overdensity(field1, pnt)
                field2 = to_overdensity(field2, pnt)
                pnt.write("calculating correlation for %s-%s"%(k1,k2))
                res = XXi(field1, field2, BOXSIZE, axis = AXIS, MAS=(MAS,MAS))

                if k1 == keylist1[0] and k2 == keylist2[0]:
                    pnt.writeTab("saving r...")
                    w.create_dataset("r",data=res.r3D)
                pnt.writeTab("saving Xi result...")
                w.create_dataset("%s-%s"%(k1,k2), data=res.xi[:,0])

                pnt.writeTab("making Xi plot...")
                plotname = PLOTS+'%s-%s/%s-%s%d_%03d'%(FILE1, FILE2, k1, k2, BOX, SNAPSHOT)
                lpt.plotCorr(res.r3D, res.xi[:,0], field1.shape[0], BOXSIZE, plotname)
    
    w.close()
    f1.close()
    f2.close()

pnt.write("2pcf.py finished")

