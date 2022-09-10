#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 15:55:52 2019

@author: Javad Vahedi & Stefan Kettemann
"""

import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
from scipy.sparse import linalg
import matplotlib.pyplot as plt
import Hamiltonian_sz
import Operator_sz
import partialtrace
import math
#import itertools
import time

#plt.rcParams['text.usetex'] = True
#plt.rcParams['text.latex.unicode'] = True
#plt.rcParams["font.family"] = "Times New Roman"

start = time.time()

# Set parameter
N  = 14#np.loadtxt('N.txt')    # number of spin
L  = 140#np.loadtxt('L.txt')   # latacice size
n  = N/L  # density
J0 = 1
hz = J0*0
sz=0
alpha = 1# np.loadtxt('alpha.txt')
n_lowest_eigenvalues = 1

f1 = open('Space_Index_L'+str(L)+'_N_'+str(N)+'_a_'+str(alpha)+'.txt', 'w')
f2 = open('Sxx_L'+str(L)+'_N_'+str(N)+'_a_'+str(alpha)+'.txt', 'w')
f3 = open('Svn_L'+str(L)+'_N_'+str(N)+'_a_'+str(alpha)+'.txt', 'w')



# Step 1: evaluate GS                                              
rows, cols, data, space_index = Hamiltonian_sz.get_hamiltonian_sparse(L, N, J0, hz, alpha,sz)
hamiltonian = sp.sparse.csr_matrix((data, (rows, cols)))
energy, groundstate = sp.sparse.linalg.eigsh(hamiltonian, 
                                             k=n_lowest_eigenvalues,
                                             which='SA',
                                             return_eigenvectors=True,
                                             maxiter=1000)

print("ground state energy (per site):", energy/N, flush=True)

del(hamiltonian)
del(energy)



groundstate = groundstate[:,0]
sxy_corrs = []
for n in range(N):
    sxy_rows, sxy_cols, sxy_data = Operator_sz.get_spinspin_sparse(L,N,0,n,hz,sz)
    sxy_matrix = sp.sparse.csr_matrix((sxy_data, (sxy_rows, sxy_cols)))
    sxy_corrs.append(np.dot(groundstate, sxy_matrix.dot(groundstate)))
    
    f2.write("%i  %i  %i  %f\n"  %(1,n+1, space_index[n] ,sxy_corrs[-1]))
    f2.flush()
f2.close()


del(sxy_matrix)


n,m,d,sxx=np.loadtxt('Sxx_L'+str(L)+'_N_'+str(N)+'_a_'+str(alpha)+'.txt',unpack=True)


#idx   = np.argsort(l)
#l = np.array(l.astype(int))[idx]
#sxx = np.array(abs(sxx))[idx]







fig = plt.figure(figsize=(4,4))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],facecolor='w')
ax.tick_params(direction='in',axis='both', color='k',labelsize=14,
               left='on', top='on', right='on', bottom='on')

ax.loglog(abs(sxx), "o-b")
x=np.arange(1,len(sxx))
ax.loglog(x,np.exp(-x),'-r')
ax.loglog(x,1/(x**alpha),'-c')

ax.axes.set_xlabel(r"$|n-m|$",fontsize=18)
ax.axes.set_ylabel(r"$\langle S^y_nS^y_m\rangle$",fontsize=18)
plt.show()


print(time.time()-start)


