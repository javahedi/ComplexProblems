
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: infer_types=True

#cimport numpy as cnp
import numpy as np
import cython
from libc.math cimport sin, cos,exp
from libc.stdio cimport printf

#@cython.boundscheck(False)
#@cython.wraparound(False)



def LorentzSmearing(double x, double x0, double sigma):
    
    '''
    Simulate the Delta function by a Lorentzian shape function
        
        \Delta(x) = \lim_{\sigma\to 0}  Lorentzian
    '''
    cdef double pi = 3.14159265359
    return 1./ pi * sigma / ((x - x0)**2 + sigma**2)


def TimeDomian( double[:] Time,  double[:] E, double[:,:] vec):
        
        cdef int dim = E.shape[0]
        cdef int nt = Time.shape[0]
        
        cdef int i =0
        cdef int j =0
        cdef int k =0
        
        
        cdef double temp = 0.0
        #cdef double complex temp
        cdef double[:] out = np.zeros(nt)
         
        #with nogil:
        for k in range(nt):
            #printf("Time  : %i\n", k)
            for i in range(dim):
                for j in range(dim):
                    out[k] += cos((E[j]-E[i])*Time[k]) * abs(vec[i,j])**2
                
        return np.array(out)



def TimeDomian2( double[:] Time,  double[:] E, double[:,:] A, double[:,:] B):
        
        cdef int dim = E.shape[0]
        cdef int nt = Time.shape[0]
        
        cdef int i =0
        cdef int j =0
        cdef int k =0
        
        
        cdef double temp = 0.0
        #cdef double complex temp
        cdef double[:] out = np.zeros(nt)
         
        #with nogil:
        for k in range(nt):
            #printf("Time  : %i\n", k)
            for i in range(dim):
                for j in range(dim):
                    out[k] += cos((E[j]-E[i])*Time[k]) * A[i,j] * B[j,i]
                
        return np.array(out)
    
    
    

    
    
    
def OmegaDomian_fixed(double[:] Omega,  double[:] E, complex[:,:] vec, double domega):
        
        cdef int dim = E.shape[0]
        cdef int nw = Omega.shape[0]
        
        cdef int i =0
        cdef int j =0
        cdef int k =0

        
        
        cdef double temp = 0.0
        #cdef double complex temp
        cdef double[:] out = np.zeros(nw)
         
        #with nogil:
        for w in range(nw):
            #printf("Omega  : %i\n", w)
            for i in range(dim):
                for j in range(dim):
                    out[w] += LorentzSmearing(Omega[w], E[j]-E[i],domega) * abs(vec[i,j])**2
                
        return np.array(out)


def OmegaDomian_log(double[:] Omega,  double[:] E, complex[:,:] vec, double[:] domega ):
        
        cdef int dim = E.shape[0]
        cdef int nw = Omega.shape[0]
        
        cdef int i =0
        cdef int j =0
        cdef int k =0

        
        
        cdef double temp = 0.0
        #cdef double complex temp
        cdef double[:] out = np.zeros(nw)
        
        #with nogil:
        for w in range(nw):
            #printf("Omega  : %i\n", w)
            for i in range(dim):
                for j in range(dim):
                    out[w]  += LorentzSmearing(Omega[w], E[j]-E[i],domega[w])     * abs(vec[i,j])**2
                
        return np.array(out)



def OmegaDomian2(double[:] Omega,  double[:] E, complex[:,:] vec, double domega, double[:] domega2 ):
        
        cdef int dim = E.shape[0]
        cdef int nw = Omega.shape[0]
        
        cdef int i =0
        cdef int j =0
        cdef int k =0

        
        
        cdef double temp = 0.0
        #cdef double complex temp
        cdef double[:] out = np.zeros(nw)
        cdef double[:] out2 = np.zeros(nw) 
        #with nogil:
        for w in range(nw):
            #printf("Omega  : %i\n", w)
            for i in range(dim):
                for j in range(dim):
                    out[w]  += LorentzSmearing(Omega[w], E[j]-E[i],domega)     * abs(vec[i,j])**2
                    out2[w] += LorentzSmearing(Omega[w], E[j]-E[i],domega2[w]) * abs(vec[i,j])**2
                
        return np.array(out),np.array(out2)