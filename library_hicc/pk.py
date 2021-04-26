#!/usr/bin/env python3

import numpy as np 
import time,sys,os
import pyfftw
import scipy.integrate as si

def frequencies(BoxSize,dims):
    kF = 2.0*np.pi/BoxSize;  middle = dims//2;  kN = middle*kF
    kmax_par = middle
    kmax_per = int(np.sqrt(middle**2 + middle**2))
    kmax     = int(np.sqrt(middle**2 + middle**2 + middle**2))
    return kF,kN,kmax_par,kmax_per,kmax

# same as frequencies but in 2D
def frequencies_2D(BoxSize,dims):
    kF = 2.0*np.pi/BoxSize;  middle = dims//2;  kN = middle*kF
    kmax_par = middle
    kmax_per = middle
    kmax     = int(np.sqrt(middle**2 + middle**2))
    return kF,kN,kmax_par,kmax_per,kmax
    
# This function finds the MAS correction index and return the array used
def MAS_function(MAS):
    MAS_index = 0;  #MAS_corr = np.ones(3,dtype=np.float64)
    if MAS=='NGP':  MAS_index = 1
    if MAS=='CIC':  MAS_index = 2
    if MAS=='TSC':  MAS_index = 3
    if MAS=='PCS':  MAS_index = 4
    return MAS_index#,MAS_corr

def MAS_correction(x, MAS_index):
    return (1.0 if (x==0.0) else pow(x/sin(x),MAS_index))

# This function performs the 3D FFT of a field in single precision
def FFT3Dr_f(a, threads)
    # align arrays
    dims  = len(a)
    a_in  = pyfftw.empty_aligned((dims,dims,dims),    dtype='float32')
    a_out = pyfftw.empty_aligned((dims,dims,dims//2+1),dtype='complex64') 

    # plan FFTW
    fftw_plan = pyfftw.FFTW(a_in, a_out, axes=(0,1,2),
                            flags=('FFTW_ESTIMATE',),
                            direction='FFTW_FORWARD', threads=threads)
                            
    # put input array into delta_r and perform FFTW
    a_in [:] = a;  fftw_plan(a_in,a_out);  return a_out

# This function performs the 3D FFT of a field in double precision
def FFT3Dr_d(a,threads):

    # align arrays
    dims  = len(a)
    a_in  = pyfftw.empty_aligned((dims,dims,dims),    dtype='float64')
    a_out = pyfftw.empty_aligned((dims,dims,dims//2+1),dtype='complex128')

    # plan FFTW
    fftw_plan = pyfftw.FFTW(a_in,a_out,axes=(0,1,2),
                            flags=('FFTW_ESTIMATE',),
                            direction='FFTW_FORWARD',threads=threads)
                            
    # put input array into delta_r and perform FFTW
    a_in [:] = a;  fftw_plan(a_in,a_out);  return a_out

# This function performs the 3D FFT of a field in single precision
def IFFT3Dr_f(a,threads):

    # align arrays
    dims  = len(a)
    a_in  = pyfftw.empty_aligned((dims,dims,dims//2+1),dtype='complex64')
    a_out = pyfftw.empty_aligned((dims,dims,dims),    dtype='float32')

    # plan FFTW
    fftw_plan = pyfftw.FFTW(a_in, a_out, axes=(0,1,2),
                            flags=('FFTW_ESTIMATE',),
                            direction='FFTW_BACKWARD', threads=threads)
                            
    # put input array into delta_r and perform FFTW
    a_in [:] = a;  fftw_plan(a_in,a_out);  return a_out

# This function performs the 3D FFT of a field in double precision
def IFFT3Dr_d(a, threads):

    # align arrays
    dims  = len(a)
    a_in  = pyfftw.empty_aligned((dims,dims,dims//2+1),dtype='complex128')
    a_out = pyfftw.empty_aligned((dims,dims,dims),    dtype='float64')

    # plan FFTW
    fftw_plan = pyfftw.FFTW(a_in,a_out,axes=(0,1,2),
                            flags=('FFTW_ESTIMATE',),
                            direction='FFTW_BACKWARD',threads=threads)
                            
    # put input array into delta_r and perform FFTW
    a_in [:] = a;  fftw_plan(a_in,a_out);  return a_out

# This function performs the 2D FFT of a field in single precision
def FFT2Dr_f(a, threads):

    # align arrays
    grid  = len(a)
    a_in  = pyfftw.empty_aligned((grid,grid),    dtype='float32')
    a_out = pyfftw.empty_aligned((grid,grid//2+1),dtype='complex64')

    # plan FFTW
    fftw_plan = pyfftw.FFTW(a_in, a_out, axes=(0,1),
                            flags=('FFTW_ESTIMATE',),
                            direction='FFTW_FORWARD', threads=threads)
                            
    # put input array into delta_r and perform FFTW
    a_in [:] = a;  fftw_plan(a_in,a_out);  return a_out

# This function performs the 2D FFT of a field in double precision
def FFT2Dr_d(a, threads):

    # align arrays
    grid  = len(a)
    a_in  = pyfftw.empty_aligned((grid,grid),    dtype='float64')
    a_out = pyfftw.empty_aligned((grid,grid//2+1),dtype='complex128')

    # plan FFTW
    fftw_plan = pyfftw.FFTW(a_in,a_out,axes=(0,1),
                            flags=('FFTW_ESTIMATE',),
                            direction='FFTW_FORWARD',threads=threads)
                            
    # put input array into delta_r and perform FFTW
    a_in [:] = a;  fftw_plan(a_in,a_out);  return a_out

# This function performs the 2D FFT of a field in single precision
def IFFT2Dr_f(a, threads):

    # align arrays
    grid  = len(a)
    a_in  = pyfftw.empty_aligned((grid,grid//2+1),dtype='complex64')
    a_out = pyfftw.empty_aligned((grid,grid),    dtype='float32')

    # plan FFTW
    fftw_plan = pyfftw.FFTW(a_in, a_out, axes=(0,1),
                            flags=('FFTW_ESTIMATE',),
                            direction='FFTW_BACKWARD', threads=threads)
                            
    # put input array into delta_r and perform FFTW
    a_in [:] = a;  fftw_plan(a_in,a_out);  return a_out

# This function performs the 2D FFT of a field in double precision
def IFFT2Dr_d(a, threads):

    # align arrays
    grid  = len(a)
    a_in  = pyfftw.empty_aligned((grid,grid//2+1),dtype='complex128')
    a_out = pyfftw.empty_aligned((grid,grid),    dtype='float64')

    # plan FFTW
    fftw_plan = pyfftw.FFTW(a_in,a_out,axes=(0,1),
                            flags=('FFTW_ESTIMATE',),
                            direction='FFTW_BACKWARD',threads=threads)
                            
    # put input array into delta_r and perform FFTW
    a_in [:] = a;  fftw_plan(a_in,a_out);  return a_out

class Pk:
    def __init__(self,delta,BoxSize,axis=2,MAS='CIC',threads=1):

        start = time.time()
        # cdef int kxx, kyy, kzz, kx, ky, kz,dims, middle, k_index, MAS_index
        # cdef int kmax_par, kmax_per, kmax, k_par, k_per, index_2D, i
        # cdef double k, delta2, prefact, mu, mu2, real, imag, kmaxper, phase
        # cdef double MAS_corr[3]
        ####### change this for double precision ######
        # cdef float MAS_factor
        # cdef np.complex64_t[:,:,::1] delta_k
        ###############################################
        # cdef np.float64_t[::1] k1D, kpar, kper, k3D, Pk1D, Pk2D, Pkphase
        # cdef np.float64_t[::1] Nmodes1D, Nmodes2D, Nmodes3D
        # cdef np.float64_t[:,::1] Pk3D 

        # find dimensions of delta: we assume is a (dims,dims,dims) array
        # determine the different frequencies and the MAS_index
        print('\nComputing power spectrum of the field...')
        dims = len(delta);  middle = dims//2
        kF,kN,kmax_par,kmax_per,kmax = frequencies(BoxSize,dims)
        MAS_index = MAS_function(MAS)

        ## compute FFT of the field (change this for double precision) ##
        delta_k = FFT3Dr_f(delta,threads)
        #################################

        # define arrays containing k1D, Pk1D and Nmodes1D. We need kmax_par+1
        # bins since modes go from 0 to kmax_par
        k1D      = np.zeros(kmax_par+1, dtype=np.float64)
        Pk1D     = np.zeros(kmax_par+1, dtype=np.float64)
        Nmodes1D = np.zeros(kmax_par+1, dtype=np.float64)

        # define arrays containing Pk2D and Nmodes2D
        Pk2D     = np.zeros((kmax_par+1)*(kmax_per+1), dtype=np.float64)
        Nmodes2D = np.zeros((kmax_par+1)*(kmax_per+1), dtype=np.float64)

        # define arrays containing k3D, Pk3D and Nmodes3D. We need kmax+1
        # bins since the mode (middle,middle, middle) has an index = kmax
        k3D      = np.zeros(kmax+1,     dtype=np.float64)
        Pk3D     = np.zeros((kmax+1,3), dtype=np.float64)
        Pkphase  = np.zeros(kmax+1,     dtype=np.float64)
        Nmodes3D = np.zeros(kmax+1,     dtype=np.float64)


        # do a loop over the independent modes.
        # compute k,k_par,k_per, mu for each mode. k's are in kF units
        start2 = time.time();  prefact = np.pi/dims
        for kxx in range(dims):
            kx = (kxx-dims if (kxx>middle) else kxx)
            MAS_corr[0] = MAS_correction(prefact*kx,MAS_index) # I'm not sure what the MAS correction index does
        
            for kyy in range(dims):
                ky = (kyy-dims if (kyy>middle) else kyy)
                MAS_corr[1] = MAS_correction(prefact*ky,MAS_index)

                for kzz in range(middle+1): #kzz=[0,1,..,middle] --> kz>0
                    kz = (kzz-dims if (kzz>middle) else kzz)
                    MAS_corr[2] = MAS_correction(prefact*kz,MAS_index)  

                    # kz=0 and kz=middle planes are special
                    if kz==0 or (kz==middle and dims%2==0):
                        if kx<0: continue
                        elif kx==0 or (kx==middle and dims%2==0):
                            if ky<0.0: continue


                    # compute |k| of the mode and its integer part
                    k       = sqrt(kx*kx + ky*ky + kz*kz)
                    k_index = int(k)

                    # compute the value of k_par and k_perp
                    if axis==0:   
                        k_par, k_per = kx, int(sqrt(ky*ky + kz*kz))
                    elif axis==1: 
                        k_par, k_per = ky, int(sqrt(kx*kx + kz*kz))
                    else:         
                        k_par, k_per = kz, int(sqrt(kx*kx + ky*ky))

                    # find the value of mu
                    if k==0:  mu = 0.0
                    else:     mu = k_par/k
                    mu2 = mu*mu

                    # take the absolute value of k_par
                    if k_par<0:  k_par = -k_par

                    # correct modes amplitude for MAS
                    MAS_factor = MAS_corr[0]*MAS_corr[1]*MAS_corr[2]
                    delta_k[kxx,kyy,kzz] = delta_k[kxx,kyy,kzz]*MAS_factor

                    # compute |delta_k|^2 of the mode
                    real = delta_k[kxx,kyy,kzz].real
                    imag = delta_k[kxx,kyy,kzz].imag
                    delta2 = real*real + imag*imag
                    phase  = atan2(real, sqrt(delta2)) 

                    # Pk1D: only consider modes with |k|<kF
                    if k<=middle:
                        k1D[k_par]      += k_par
                        Pk1D[k_par]     += delta2
                        Nmodes1D[k_par] += 1.0

                    # Pk2D: P(k_per,k_par)
                    # index_2D goes from 0 to (kmax_par+1)*(kmax_per+1)-1
                    index_2D = (kmax_par+1)*k_per + k_par
                    Pk2D[index_2D]     += delta2
                    Nmodes2D[index_2D] += 1.0

                    # Pk3D.
                    k3D[k_index]      += k
                    Pk3D[k_index,0]   += delta2
                    Pk3D[k_index,1]   += (delta2*(3.0*mu2-1.0)/2.0)
                    Pk3D[k_index,2]   += (delta2*(35.0*mu2*mu2 - 30.0*mu2 + 3.0)/8.0)
                    Pkphase[k_index]  += (phase*phase)
                    Nmodes3D[k_index] += 1.0
        print('Time to complete loop = %.2f'%(time.time()-start2))

        # Pk1D. Discard DC mode bin and give units
        # the perpendicular modes sample an area equal to pi*kmax_per^2
        # we assume that each mode has an area equal to pi*kmax_per^2/Nmodes
        k1D  = k1D[1:];  Nmodes1D = Nmodes1D[1:];  Pk1D = Pk1D[1:]
        for i in range(len(k1D)):
            Pk1D[i] = Pk1D[i]*(BoxSize/dims**2)**3 #give units
            k1D[i]  = (k1D[i]/Nmodes1D[i])*kF      #give units
            kmaxper = sqrt(kN**2 - k1D[i]**2)
            Pk1D[i] = Pk1D[i]*(np.pi*kmaxper**2/Nmodes1D[i])/(2.0*np.pi)**2
        self.k1D = np.asarray(k1D);  self.Pk1D = np.asarray(Pk1D)
        self.Nmodes1D = np.asarray(Nmodes1D  )

        # Pk2D. Keep DC mode bin, give units to Pk2D and find kpar & kper
        kpar = np.zeros((kmax_par+1)*(kmax_per+1), dtype=np.float64)
        kper = np.zeros((kmax_par+1)*(kmax_per+1), dtype=np.float64)
        for k_par in range(kmax_par+1):
            for k_per in range(kmax_per+1):
                index_2D = (kmax_par+1)*k_per + k_par
                kpar[index_2D] = 0.5*(k_par + k_par+1)*kF
                kper[index_2D] = 0.5*(k_per + k_per+1)*kF
        for i in range(len(kpar)):
            Pk2D[i] = Pk2D[i]*(BoxSize/dims**2)**3/Nmodes2D[i]
        self.kpar = np.asarray(kpar);  self.kper = np.asarray(kper)
        self.Pk2D = np.asarray(Pk2D);  self.Nmodes2D = np.asarray(Nmodes2D)

        # Pk3D. Check modes, discard DC mode bin and give units
        # we need to multiply the multipoles by (2*ell + 1)
        check_number_modes(Nmodes3D,dims)
        k3D  = k3D[1:];  Nmodes3D = Nmodes3D[1:];  Pk3D = Pk3D[1:,:]
        Pkphase = Pkphase[1:]
        for i in range(len(k3D)):
            k3D[i]     = (k3D[i]/Nmodes3D[i])*kF
            Pk3D[i,0]  = (Pk3D[i,0]/Nmodes3D[i])*(BoxSize/dims**2)**3
            Pk3D[i,1]  = (Pk3D[i,1]*5.0/Nmodes3D[i])*(BoxSize/dims**2)**3
            Pk3D[i,2]  = (Pk3D[i,2]*9.0/Nmodes3D[i])*(BoxSize/dims**2)**3
            Pkphase[i] = (Pkphase[i]/Nmodes3D[i])*(BoxSize/dims**2)**3
        self.k3D = np.asarray(k3D);  self.Nmodes3D = np.asarray(Nmodes3D)
        self.Pk = np.asarray(Pk3D);  self.Pkphase = Pkphase

        print('Time taken = %.2f seconds'%(time.time()-start))
        return

class XPk:
    def __init__(self,delta,BoxSize,axis=2,MAS=None,threads=1):

        start = time.time()
        # cdef int dims, middle, fields, Xfields, num_unique_MAS, i, j
        # cdef int index, index_z, index_2D, index_X
        # cdef int kxx, kyy, kzz, kx, ky, kz, k_index#, begin, end
        # cdef int kmax_par, kmax_per, kmax, k_par, k_per
        # cdef double k, prefact, mu, mu2, val1, val2
        # cdef double delta2, delta2_X, fact
        ####### change this for double precision ######
        # cdef float MAS_factor
        # cdef np.ndarray[np.complex64_t,ndim=3] delta_k
        #cdef np.complex64_t[:,:,::1] delta_k
        ###############################################
        # cdef np.int32_t[::1] MAS_index, unique_MAS_id
        # cdef np.float64_t[::1] real_part, imag_part, k1D, k3D
        # cdef np.float64_t[::1] kpar, kper, Nmodes1D, Nmodes2D, Nmodes3D
        # cdef np.float64_t[:,::1] MAS_corr, Pk1D, Pk2D, PkX1D, PkX2D
        # cdef np.float64_t[:,:,::1] Pk3D, PkX3D

        print('\nComputing power spectra of the fields...')

        # find the number and dimensions of the density fields
        # we assume the density fields are (dims,dims,dims) arrays
        dims = len(delta[0]);  middle = dims//2;  fields = len(delta)
        Xfields = fields*(fields-1)//2  #number of independent cross-P(k)

        # check that the dimensions of all fields are the same
        for i in range(1,fields):
            if len(delta[i])!=dims:
                print('Fields have different grid sizes!!!'); sys.exit()

        # find the different relevant frequencies
        kF,kN,kmax_par,kmax_per,kmax = frequencies(BoxSize,dims)

        # find the independent MAS and the arrays relating both.
        # if MAS = ['CIC','PCS','CIC','CIC'] ==> unique_MAS = ['CIC','PCS']
        # num_unique_MAS = 2 : unique_MAS_id = [0,1,0,0]
        unique_MAS     = np.array(list(set(MAS))) #array with independent MAS
        num_unique_MAS = len(unique_MAS)          #number of independent MAS
        unique_MAS_id  = np.empty(fields,dtype=np.int32) 
        for i in range(fields):
            unique_MAS_id[i] = np.where(MAS[i]==unique_MAS)[0][0]

        # define and fill the MAS_corr and MAS_index arrays
        MAS_corr  = np.ones((num_unique_MAS,3), dtype=np.float64)
        MAS_index = np.zeros(num_unique_MAS,    dtype=np.int32)
        for i in range(num_unique_MAS):
            MAS_index[i] = MAS_function(unique_MAS[i])

        # define the real_part and imag_part arrays
        real_part = np.zeros(fields,dtype=np.float64)
        imag_part = np.zeros(fields,dtype=np.float64)

        ## compute FFT of the field (change this for double precision) ##
        # to try to have the elements of the different fields as close as 
        # possible we stack along the z-direction (major-row)
        delta_k = np.empty((dims,dims,(middle+1)*fields),dtype=np.complex64)
        for i in range(fields):
            begin = i*(middle+1);  end = (i+1)*(middle+1)
            delta_k[:,:,begin:end] = FFT3Dr_f(delta[i],threads)
        print('Time FFTS = %.2f'%(time.time()-start))
        #################################

        # define arrays having k1D, Pk1D, PkX1D & Nmodes1D. We need kmax_par+1
        # bins since modes go from 0 to kmax_par. Is better if we define the
        # arrays as (kmax_par+1,fields) rather than (fields,kmax_par+1) since
        # in memory arrays numpy arrays are row-major
        k1D      = np.zeros(kmax_par+1,           dtype=np.float64)
        Pk1D     = np.zeros((kmax_par+1,fields),  dtype=np.float64)
        PkX1D    = np.zeros((kmax_par+1,Xfields), dtype=np.float64)
        Nmodes1D = np.zeros(kmax_par+1,           dtype=np.float64)

        # define arrays containing Pk2D and Nmodes2D. We define the arrays
        # in this way to have them as close as possible in row-major
        Pk2D     = np.zeros(((kmax_par+1)*(kmax_per+1),fields), 
                            dtype=np.float64)
        PkX2D    = np.zeros(((kmax_par+1)*(kmax_per+1),Xfields),
                            dtype=np.float64)
        Nmodes2D = np.zeros((kmax_par+1)*(kmax_per+1), 
                            dtype=np.float64)

        # define arrays containing k3D, Pk3D,PkX3D & Nmodes3D. We need kmax+1
        # bins since the mode (middle,middle, middle) has an index = kmax.
        # We define the arrays in this way to benefit of row-major
        k3D      = np.zeros(kmax+1,             dtype=np.float64)
        Pk3D     = np.zeros((kmax+1,3,fields),  dtype=np.float64)
        PkX3D    = np.zeros((kmax+1,3,Xfields), dtype=np.float64)
        Nmodes3D = np.zeros(kmax+1,             dtype=np.float64)

        # do a loop over the independent modes.
        # compute k,k_par,k_per, mu for each mode. k's are in kF units
        start2 = time.time();  prefact = np.pi/dims
        for kxx in range(dims):
            kx = (kxx-dims if (kxx>middle) else kxx)
            for i in range(num_unique_MAS):
                MAS_corr[i,0] = MAS_correction(prefact*kx,MAS_index[i])

            for kyy in range(dims):
                ky = (kyy-dims if (kyy>middle) else kyy)
                for i in range(num_unique_MAS):
                    MAS_corr[i,1] = MAS_correction(prefact*ky,MAS_index[i])

                for kzz in range(middle+1): #kzz=[0,1,..,middle] --> kz>0
                    kz = (kzz-dims if (kzz>middle) else kzz)
                    for i in range(num_unique_MAS):
                        MAS_corr[i,2] = MAS_correction(prefact*kz,MAS_index[i])

                    # kz=0 and kz=middle planes are special
                    if kz==0 or (kz==middle and dims%2==0):
                        if kx<0: continue
                        elif kx==0 or (kx==middle and dims%2==0):
                            if ky<0.0: continue

                    ###### k, k_index, k_par,k_per, mu ######
                    # compute |k| of the mode and its integer part
                    k       = sqrt(kx*kx + ky*ky + kz*kz)
                    k_index = int(k)

                    # compute the value of k_par and k_perp
                    if axis==0:   
                        k_par, k_per = kx, int(sqrt(ky*ky + kz*kz))
                    elif axis==1: 
                        k_par, k_per = ky, int(sqrt(kx*kx + kz*kz))
                    else:         
                        k_par, k_per = kz, int(sqrt(kx*kx + ky*ky))

                    # find the value of mu
                    if k==0:  mu = 0.0
                    else:     mu = k_par/k
                    mu2 = mu*mu  
                    val1 = (3.0*mu2-1.0)/2.0
                    val2 = (35.0*mu2*mu2 - 30.0*mu2 + 3.0)/8.0

                    # take the absolute value of k_par
                    if k_par<0:  k_par = -k_par
                    #########################################

                    ####### fill the general arrays #########
                    # Pk1D(k)
                    if k<=middle:
                        k1D[k_par]      += k_par
                        Nmodes1D[k_par] += 1.0
                    
                    # Pk2D: index_2D goes from 0 to (kmax_par+1)*(kmax_per+1)-1
                    index_2D = (kmax_par+1)*k_per + k_par
                    Nmodes2D[index_2D] += 1.0

                    # Pk3D
                    k3D[k_index]      += k
                    Nmodes3D[k_index] += 1.0
                    #########################################

                    #### correct modes amplitude for MAS ####
                    for i in range(fields):
                        index = unique_MAS_id[i]
                        MAS_factor = MAS_corr[index,0]*\
                                        MAS_corr[index,1]*\
                                        MAS_corr[index,2]
                        index_z = i*(middle+1) + kzz
                        delta_k[kxx,kyy,index_z] = delta_k[kxx,kyy,index_z]*\
                                                    MAS_factor
                        real_part[i] = delta_k[kxx,kyy,index_z].real
                        imag_part[i] = delta_k[kxx,kyy,index_z].imag

                        ########## compute auto-P(k) ########
                        delta2 = real_part[i]*real_part[i] +\
                                    imag_part[i]*imag_part[i]

                        # Pk1D: only consider modes with |k|<kF
                        if k<=middle:
                            Pk1D[k_par,i] += delta2

                        # Pk2D: P(k_per,k_par)
                        Pk2D[index_2D,i] += delta2

                        # Pk3D
                        Pk3D[k_index,0,i] += (delta2)
                        Pk3D[k_index,1,i] += (delta2*val1)
                        Pk3D[k_index,2,i] += (delta2*val2)
                    #########################################

                    ####### compute XPk for each pair #######
                    index_X  = 0
                    for i in range(fields):
                        for j in range(i+1,fields):
                            delta2_X = real_part[i]*real_part[j] +\
                                        imag_part[i]*imag_part[j]            

                            # Pk1D: only consider modes with |k|<kF
                            if k<=middle:
                                PkX1D[k_par,index_X] += delta2_X

                            # Pk2D: P(k_per,k_par)
                            PkX2D[index_2D,index_X] += delta2_X
                            
                            # Pk3D
                            PkX3D[k_index,0,index_X] += delta2_X
                            PkX3D[k_index,1,index_X] += (delta2_X*val1)
                            PkX3D[k_index,2,index_X] += (delta2_X*val2)

                            index_X += 1
                    #########################################

        print('Time loop = %.2f'%(time.time()-start2))
        fact = (BoxSize/dims**2)**3

        # Pk1D. Discard DC mode bin and give units
        # the perpendicular modes sample an area equal to pi*kmax_per^2
        # we assume that each mode has an area equal to pi*kmax_per^2/Nmodes
        k1D  = k1D[1:];  Nmodes1D = Nmodes1D[1:]
        Pk1D = Pk1D[1:,:];  PkX1D = PkX1D[1:,:]
        for i in range(len(k1D)):
            k1D[i]  = (k1D[i]/Nmodes1D[i])*kF  #give units
            kmaxper = sqrt(kN**2 - k1D[i]**2)

            for j in range(fields):
                Pk1D[i,j] = Pk1D[i,j]*fact #give units
                Pk1D[i,j] = Pk1D[i,j]*(np.pi*kmaxper**2/Nmodes1D[i])/(2.0*np.pi)**2
            for j in range(Xfields):
                PkX1D[i,j] = PkX1D[i,j]*fact #give units
                PkX1D[i,j] = PkX1D[i,j]*(np.pi*kmaxper**2/Nmodes1D[i])/(2.0*np.pi)**2
        self.k1D = np.asarray(k1D);    self.Nmodes1D = np.asarray(Nmodes1D)  
        self.Pk1D = np.asarray(Pk1D);  self.PkX1D = np.asarray(PkX1D)

        # Pk2D. Keep DC mode bin, give units to Pk2D and find kpar & kper
        kpar = np.zeros((kmax_par+1)*(kmax_per+1), dtype=np.float64)
        kper = np.zeros((kmax_par+1)*(kmax_per+1), dtype=np.float64)
        for k_par in range(kmax_par+1):
            for k_per in range(kmax_per+1):
                index_2D = (kmax_par+1)*k_per + k_par
                kpar[index_2D] = 0.5*(k_par + k_par+1)*kF
                kper[index_2D] = 0.5*(k_per + k_per+1)*kF
        for i in range(len(kpar)):
            for j in range(fields):
                Pk2D[i,j] = Pk2D[i,j]*fact/Nmodes2D[i]
            for j in range(Xfields):
                PkX2D[i,j] = PkX2D[i,j]*fact/Nmodes2D[i]
        self.kpar = np.asarray(kpar);  self.kper = np.asarray(kper)
        self.Nmodes2D = np.asarray(Nmodes2D)
        self.Pk2D = np.asarray(Pk2D);  self.PkX2D = np.asarray(PkX2D)

        # Pk3D. Check modes, discard DC mode bin and give units
        # we need to multiply the multipoles by (2*ell + 1)
        check_number_modes(Nmodes3D,dims)
        k3D  = k3D[1:];  Nmodes3D = Nmodes3D[1:];  
        Pk3D = Pk3D[1:,:,:];  PkX3D = PkX3D[1:,:,:]
        for i in range(len(k3D)):
            k3D[i] = (k3D[i]/Nmodes3D[i])*kF

            for j in range(fields):
                Pk3D[i,0,j] = (Pk3D[i,0,j]/Nmodes3D[i])*fact
                Pk3D[i,1,j] = (Pk3D[i,1,j]*5.0/Nmodes3D[i])*fact
                Pk3D[i,2,j] = (Pk3D[i,2,j]*9.0/Nmodes3D[i])*fact

            for j in range(Xfields):
                PkX3D[i,0,j] = (PkX3D[i,0,j]/Nmodes3D[i])*fact
                PkX3D[i,1,j] = (PkX3D[i,1,j]*5.0/Nmodes3D[i])*fact
                PkX3D[i,2,j] = (PkX3D[i,2,j]*9.0/Nmodes3D[i])*fact

        self.k3D = np.asarray(k3D);  self.Nmodes3D = np.asarray(Nmodes3D)
        self.Pk = np.asarray(Pk3D);  self.XPk = np.asarray(PkX3D)

        print('Time taken = %.2f seconds'%(time.time()-start))
        return