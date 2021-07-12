#!/usr/bin/env python3
"""
Calculates the power spectrum of the fields given by a text file.
The name of the text file is given as command-line input, .
"""

# import statements
import sys
import numpy as np
import h5py as hp
import time
import illustris_python as il
from Pk_library import Pk
from Pk_library import XPk

# reading command line inputs
start = time.time()
AUTO_OR_XPK = sys.argv[1]
SNAPSHOT = int(sys.argv[2])
BOX = int(sys.argv[3])
DIM = int(sys.argv[4])
AXIS = int(sys.argv[5])
FILE1 = sys.argv[6]
if AUTO_OR_XPK == "cross":
    FILE2 = sys.argv[7]
    IS_XPK = True
    SAME_FILE = FILE1 == FILE2
elif AUTO_OR_XPK == "auto":
    IS_XPK = False
    SMALL_GRID = 0
else:
    raise ValueError("incorrect input: first arg should be auto or cross")

print("read command-line inputs")
# defining basic paths
HOME = '/lustre/cosinga/final_fields/'
TNG = '/lustre/cosinga/tng%d'%BOX

# getting author-defined constants
MAS = 'CIC'

# memory issues, create log file to write to
if IS_XPK:
    log = open("/lustre/cosinga/hicc/pk/%s-%spk_log.log"%(FILE1,FILE2),'w')
else:
    log = open("/lustre/cosinga/hicc/pk/%spk_log.log"%(FILE1),'w')

# getting simulation defined constants
head = il.groupcat.loadHeader(TNG, SNAPSHOT)
SCALE = head['Time'] # scale factor
BOXSIZE = head['BoxSize']/1e3 * SCALE #Mpc/h - convert from comoving using a
del head
print("The boxsize is %.3f"%BOXSIZE)
print("got simulation-defined constants")
def to_overdensity(field):
    if SMALL_GRID == 1:
        field = field[::2,:,:]+field[1::2,:,:]
        field = field[:, ::2, :] + field[:, ::2, :]
        field = field[:, :, ::2] + field[:, :, ::2]
    else:
        print("not halving the size of the grid")
    field = field/BOXSIZE**3
    field = field/np.mean(field).astype(np.float32)
    field=field - 1
    if not field.dtype == np.float32:
        print("turning into float32s")
        field.astype(np.float32)
    
    return field
    

if IS_XPK:
    print("calculating the cross-power for %s, %s"%(FILE1, FILE2))

    log.write("calculating the cross-power for %s, %s\n"%(FILE1, FILE2))
    # output file
    if DIM==0: # create both 1D and 2D
        w1 = hp.File(HOME+'pk/%s-%s%d_%03d.1Dxpk.hdf5'%(FILE1,FILE2,BOX,SNAPSHOT),'w')
        w2 = hp.File(HOME+'pk/%s-%s%d_%03d.2Dxpk.hdf5'%(FILE1,FILE2,BOX,SNAPSHOT),'w')
    else:
        w = hp.File(HOME+'pk/%s-%s%d_%03d.%dDxpk.hdf5'%(FILE1,FILE2,BOX,SNAPSHOT,DIM),'w')

    # input files, loop over all of the keys in each
    f1 = hp.File(HOME+FILE1+'%d_%03d.final.hdf5'%(BOX, SNAPSHOT),'r')
    key1list = list(f1.keys())
    f2 = hp.File(HOME+FILE2+'%d_%03d.final.hdf5'%(BOX, SNAPSHOT),'r')
    key2list = list(f2.keys())
    log.write('opened both files, size of each are %.3e and %.3e.\n starting calculations... \n'%(sys.getsizeof(f1),sys.getsizeof(f2)))
    for key1 in key1list:
        for key2 in key2list:
            log.write("now calculating xpk for %s-%s...\n"%(key1,key2))
            # get the fields
            field1 = f1[key1][:]
            log.write("committed the first field, size=%.3e\n"%sys.getsizeof(field1))
            field2 = f2[key2][:]
            log.write("committed the second field, size=%.3e\n"%sys.getsizeof(field2))
            if SAME_FILE and key1 == key2:
                print("skipping xpk calculation for %s, %s; just auto power"%(key1, key2))
                log.write("skipping xpk calculation for %s, %s; just auto power\n"%(key1,key2))
            elif not (len(field1.shape) == 3 or len(field2.shape)==3):
                print("skipping calculation for %s, %s; grid is not correct shape"%(key1, key2))
            else:
                print("starting xpk calculation for %s, %s"%(key1, key2))
                
                # convert them to overdensities
                field1=to_overdensity(field1)
                field2=to_overdensity(field2)
                # compute the xpk
                res = XPk([field1,field2], BOXSIZE, axis=AXIS, MAS=[MAS, MAS])
                print("just got the xpk result, size = %.3e\n"%sys.getsizeof(res))
                # if this is the first calculation, save the wavenumbers
                if key1 == key1list[0] and key2 == key2list[0]:
                    # if in 1D, just need 1 k, otherwise need kper and kpar
                    if DIM == 1:
                        w.create_dataset("k",data=res.k3D)
                    elif DIM == 2:
                        w.create_dataset("kper", data=res.kper)
                        w.create_dataset("kpar", data=res.kpar)
                    elif DIM == 0:
                        w1.create_dataset("k",data=res.k3D)
                        w2.create_dataset("kper", data=res.kper)
                        w2.create_dataset("kpar", data=res.kpar)            
                # save xpk result
                log.write("now saving the xpk result...\n")
                if DIM == 1:
                    w.create_dataset("%s-%s"%(key1, key2),data=res.XPk[:,0,0])
                elif DIM == 2:
                    w.create_dataset("%s-%s"%(key1, key2), data=res.PkX2D[:,0])
                elif DIM == 0:
                    w1.create_dataset("%s-%s"%(key1, key2),data= res.XPk[:,0,0])
                    w2.create_dataset("%s-%s"%(key1, key2),data=res.PkX2D[:,0])
                    # 2D only has a field index, no "ell" index
    f1.close()
    f2.close()
    
    if DIM==0:
        
        w1.close()
        w2.close()
    else:
        w.close()
else:# auto power spectrum
    print("starting procedure for the auto power for %s"%FILE1)
    # output file
    if DIM==0: # create both 1D and 2D
        w1 = hp.File(HOME+'pk/%s%d_%03d.1Dpk.hdf5'%(FILE1,BOX,SNAPSHOT),'w')
        w2 = hp.File(HOME+'pk/%s%d_%03d.2Dpk.hdf5'%(FILE1,BOX,SNAPSHOT),'w')
    else:
        w = hp.File(HOME+'pk/%s%d_%03d.%dDpk.hdf5'%(FILE1,BOX,SNAPSHOT,DIM),'w')

    # input data
    f1 = hp.File(HOME+FILE1+'%d_%03d.final.hdf5'%(BOX, SNAPSHOT),'r')
    keylist = list(f1.keys())
    
    # iterate over each k, save its auto power to file
    for k in keylist:

        field1 = f1[k][:]
        if not len(field1.shape) == 3:
            print("skipping pk calculation for %s; not the correct shape."%k)
        else:
            print("calculating pk for %s"%k)
            field1=to_overdensity(field1)
            res = Pk(field1, BOXSIZE, axis=AXIS, MAS=MAS)

            # if first calculation, save the wavenumbers
            if k == keylist[0]:
                print("saving the wavenumbers...")
                # if in 1D, just need 1 k, otherwise need kper and kpar
                if DIM == 1:
                    w.create_dataset("k",data=res.k3D)
                elif DIM == 2:
                    w.create_dataset("kper", data=res.kper)
                    w.create_dataset("kpar", data=res.kpar)
                elif DIM == 0:
                    w1.create_dataset("k",data=res.k3D)
                    w2.create_dataset("kper", data=res.kper)
                    w2.create_dataset("kpar", data=res.kpar)
            
            # save pk result
            print("saving the pk result...")
            if DIM == 1:
                w.create_dataset(k,data= res.Pk[:,0])
            elif DIM == 2:
                w.create_dataset(k,data=res.Pk2D[:])
            elif DIM == 0:
                w1.create_dataset(k,data= res.Pk[:,0])
                w2.create_dataset(k,data=res.Pk2D[:])
    f1.close()
    if DIM==0:
        w1.close()
        w2.close()
    else:
        runtime = time.time()-start
        w.create_dataset("runtime", data=runtime)
        w.close()
    print("\nfinished power spectrum calculation.")
