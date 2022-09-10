#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 15:55:52 2019
@author: Javad Vahedi & Stefan Kettemann

"""

import sys,os
#os.environ['OMP_NUM_THREADS']='2' # set number of OpenMP threads to run in parallel
#os.environ['MKL_NUM_THREADS']='2' # set number of MKL threads to run in parallel



import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
from scipy.sparse import linalg
import Hamiltonian
import math
import Operator
import time
import matplotlib.pylab as plt
import random

# Set parameter
N  = 12#int(np.loadtxt('N.txt'))    # number of spin
L  = 120#int(np.loadtxt('L.txt'))    # lattice size

J0 = 1
hz = J0*0
alpha =  1.0#np.loadtxt('alpha.txt')
n_lowest_eigenvalues = 1

nrandom  = 1000#int(np.loadtxt('random.txt'))    # number of random

for w in range(100,nrandom):

    start = time.time()
    
    f1 = open('Sxx_'+str(w)+'.txt', 'w')
    #f2 = open('Svn_SVD_'+str(w)+'.txt','w')
    f3 = open('Svn_RDM_site_'+str(w)+'.txt', 'w')
    f4 = open('Svn_RDM_length_'+str(w)+'.txt', 'w')
 

    site_index =random.sample(range(L), N)
    site_index.sort()
    
    # Step 1: evaluate GS
    rows, cols, data = Hamiltonian.get_hamiltonian_sparse(L, N, J0, hz, site_index,alpha,bc='obc')
    hamiltonian = sp.sparse.csr_matrix((data, (rows, cols)))

    #clear memory
    del rows, cols, data

    energy, groundstate = sp.sparse.linalg.eigsh(hamiltonian,
                                                 k=n_lowest_eigenvalues,
                                                 which='SA',
                                                 return_eigenvectors=True,
                                                 maxiter=1E4)

    #clear memory
    del hamiltonian
    groundstate = groundstate[:,0]

    #%%s
    ###################################
    #  Sxx correlation
    ####################################
    sxy_corrs = []
    for m in range(N):
        for n in range(m+1,N):
                d = site_index[n] - site_index[m]
                sxy_rows, sxy_cols, sxy_data = Operator.get_spinspin_sparse(L,N,m,n,hz,False)
                sxy_matrix = sp.sparse.csr_matrix((sxy_data, (sxy_rows, sxy_cols)))
                sxy_corrs.append(np.dot(groundstate, sxy_matrix.dot(groundstate)))
                f1.write("%i  %i  %i  %f   %f\n"  %(m,n, d, sxy_corrs[-1],abs(sxy_corrs[-1])))
                f1.flush()
    f1.close()

    del sxy_matrix, sxy_rows, sxy_cols, sxy_data

    '''
    #%%
    ###################################
    #  Von Numann entropy with
    #  Singula value decomposition (SVD)
    ####################################
    en = []
    for i in range(1,N):
        #start = time.time()
        GS = np.reshape(groundstate,(int(2**i),int(2**(N-i))))
        s = sp.sparse.linalg.svds(GS,
                              min(GS.shape)-1,
                              which='LM',
                              return_singular_vectors=False)
        
        p=abs(s)**2
        p[p <= 0] = 1e-20
        en.append(np.sum(-p*np.log(p)))
        #print('Across bond b=%i, time=%f,  SvN =%f'%(i,time.time()-start, en[-1]), flush=True)
        
        f2.write("%i  %f\n"  %(i, en[-1]))
        f2.flush()

    f2.close()
    '''



    #%%
    ###################################
    #  Von Numann entropy with
    #  Reduced Density Matrix (RDM)
    ####################################

    SvnL = []
    SvnR = []
    for i in range(1,N//2+1):
          
        da = 2**(N-i)
        db = 2**i
        
        rho = np.reshape(groundstate,(da,db)).conj().T @ np.reshape(groundstate,(da,db))
        eigs = sp.linalg.eigh(rho,eigvals_only=True)
        del(rho)

        eigs[eigs <= 0] = 1e-20
        svn  = -np.sum(np.real(eigs*np.log(eigs)))
        if math.isnan(svn):
            SvnL.append(0)
        else:
            SvnL.append(-np.sum(np.real(eigs*np.log(eigs))))
        #print('Across bond b=%i, time=%f,  SvN =%f'%(i,time.time()-start, SvnL[-1]), flush=True)
        

        rho = np.reshape(groundstate,(db,da)) @ np.reshape(groundstate,(db,da)).conj().T
        eigs = sp.linalg.eigh(rho,eigvals_only=True)
        del(rho)

        eigs[eigs <= 0] = 1e-20
        svn  = -np.sum(np.real(eigs*np.log(eigs)))
        if math.isnan(svn):
            SvnR.append(0)
        else:
            SvnR.append(svn)
        #print('Across bond b=%i, time=%f,  SvN =%f'%(i,time.time()-start, SvnR[-1]), flush=True)
        #print('                         ')

    Svn = SvnL+SvnR[::-1][1::]

    for i in Svn:
        f3.write('%f\n'%(i))
    f3.close()

    ##########################
    ##   SvN  based on length
    ##########################
    SvnTotal = [0]
    for l in range(L):
        if l in site_index:
            indx = np.argwhere(np.asanyarray(site_index)==l)[0][0]
            SvnTotal.append(Svn[indx-1])
        else:
            SvnTotal.append(SvnTotal[-1])
        
        if l>site_index[-1]:
           SvnTotal.append(0)

        f4.write('%i   %f\n'%(l,SvnTotal[-1]))
        f4.flush()
    f4.close()
    
    final = time.time()-start
    output = f" loop number is {w} and time needs is {final} second."
    print(output, flush=True)
