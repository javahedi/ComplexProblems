#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 15:55:52 2019

@author: Javad Vahedi & Stefan Kettemann


python -m numpy.f2py3 -c partial_trace.f90 partial_trace

or 

f2py3 -c -m  partial_trace partial_trace.f90 


import partial_trace
print(partial_trace.__doc__)
"""



import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
from scipy.sparse import linalg
import matplotlib.pyplot as plt
import Hamiltonian
import Operator
#import partialtrace
import math
#import ptrace
import ptrace_saeed
import ptrace_SAEED

#import partial_trace
#import Concurence
#import entanglement
import time

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True
plt.rcParams['font.family'] = 'Times New Roman'




start = time.time()

# Set parameter
N  = 10      #np.loadtxt('N.txt')    # number of spin
L  = 100     #np.loadtxt('L.txt')    # lattice size
n  = N/L     # density
J0 = 1
hz = J0*0
alpha = 1    # np.loadtxt('alpha.txt')
eta = 9
n_lowest_eigenvalues = 1

#f1 = open('Space_Index_L'+str(L)+'_N_'+str(N)+'_a_'+str(alpha)+'.txt', 'w')
#f2 = open('Sxx_L'+str(L)+'_N_'+str(N)+'_a_'+str(alpha)+'.txt', 'w')
#f3 = open('Svn_L'+str(L)+'_N_'+str(N)+'_a_'+str(alpha)+'.txt', 'w')



# Step 1: evaluate GS 
rows, cols, data  = Hamiltonian.get_hamiltonian_sparse(L, N, J0, hz, alpha,eta)
hamiltonian = sp.sparse.csr_matrix((data, (rows, cols)))
energy, groundstate = sp.sparse.linalg.eigsh(hamiltonian, 
                                             k=n_lowest_eigenvalues,
                                             which='SA',
                                             return_eigenvectors=True,            
                                             maxiter=1000)

print("ground state energy (per site):", energy[0]/N, flush=True)
groundstate = groundstate[:,0]
groundstate[abs(groundstate) < 1e-8] = 0

#Psi0=sp.sparse.csr_matrix(groundstate)

del(hamiltonian)
del(energy)

'''
sxy_corrs = []
for n in range(1,N):
    sxy_rows, sxy_cols, sxy_data = Operator.get_spinspin_sparse(L,N,0,n,hz,False)
    sxy_matrix = sp.sparse.csr_matrix((sxy_data, (sxy_rows, sxy_cols)))
    sxy_corrs.append(np.dot(groundstate, sxy_matrix.dot(groundstate)))
    f2.write("%i  %i  %i  %f\n"  %(1,n+1, space_index[n] ,sxy_corrs[-1]))
    f2.flush()
    
    print(sxy_corrs[-1])
f2.close()


del(sxy_matrix)
'''

 
    


####################################
#  Von Numann entropy
####################################
nss = N         # Number of subsystems
d   = 2**N      # Total dimension
di  = [2]*N     # Vector specifying the dimensions of the sub-systems
Svn = []
for i in range(1,N):
    
   
    
    
    db = 2**i
    da = 2**(N-i)
    
    #rho = ptrace_saeed.get_trace(groundstate, da, db)
    #eigs = np.linalg.eigvalsh(rho)
    

    rho = ptrace_SAEED.get_trace(groundstate, da, db)
    eigs = sp.linalg.eigh(rho.toarray(),eigvals_only=True)                                  
    
    eigs[eigs <= 0] = 1e-20
    svn  = -np.sum(np.real(eigs*np.log(eigs)))
    if math.isnan(svn):
        Svn.append(0)
    else:
        Svn.append(-np.sum(np.real(eigs*np.log(eigs))))
    print('Across bond b=%i, SvN =%f'%(i,Svn[-1]))
    
    
    
plt.plot(Svn,'o-r')    

    
   

"""
for i in range(1,N+1):
    a=[0]*i+[1]*(N-i)
    b=a[::-1]
    print(a+[1]+b)
"""
    
    
    






