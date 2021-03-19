#!/usr/bin/env python3
"""
Moved Paco's MAS function to here as a temporary workaround.
"""
import numpy as np

def CICW(pos, number, BoxSize, mass):
    ptls = pos.shape[0]; coord = pos.shape[1]; dims = number.shape[0]
    inv_cell_size = dims/BoxSize
    
    index_d = np.zeros(3)
    index_u = np.zeros(3)
    d = np.zeros(3)
    u = np.zeros(3)

    for i in range(ptls):
        for axis in range(coord):
            dist = pos[i,axis] * inv_cell_size
            u[axis] = dist - int(dist)
            d[axis] = 1 - u[axis]
            index_d[axis] = (int(dist))%dims
            index_u[axis] = index_d[axis] + 1
            index_u[axis] = index_u[axis]%dims #seems this is faster
        number[index_d[0],index_d[1],index_u[2]] += d[0]*d[1]*u[2]
        number[index_d[0],index_u[1],index_d[2]] += d[0]*u[1]*d[2]
        number[index_d[0],index_u[1],index_u[2]] += d[0]*u[1]*u[2]
        number[index_u[0],index_d[1],index_d[2]] += u[0]*d[1]*d[2]
        number[index_u[0],index_d[1],index_u[2]] += u[0]*d[1]*u[2]
        number[index_u[0],index_u[1],index_d[2]] += u[0]*u[1]*d[2]
        number[index_u[0],index_u[1],index_u[2]] += u[0]*u[1]*u[2]