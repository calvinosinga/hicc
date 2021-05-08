import numpy as np
import time
def pos_redshift_space(pos, vel, Boxsize, Hubble, redshift, axis):
    start = time.time()
    print("shifting to redshift space...")
    factor = (1+redshift)/Hubble
    pos[:,axis] = pos[:,axis] + vel[:,axis]*factor
    pos[:,axis] = np.where((pos[:,axis]>Boxsize) | (pos[:,axis]<0), (pos[:,axis]+Boxsize)%Boxsize, pos[:,axis])
    # if redshift space pushes some particles beyond the boundary of the box, use periodic assumption
    print("took %.3e seconds"%(time.time()-start))
    return pos
    