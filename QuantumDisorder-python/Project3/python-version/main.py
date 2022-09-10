#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  7 14:57:07 2019

@author: javadvahedi
"""

import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
from scipy.sparse import linalg
import matplotlib.pyplot as plt
import Hamiltonian
#import Operator
import time


import matplotlib as mpl
mpl.rcParams['text.usetex']        = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.family']        = 'sans-serif'
mpl.rcParams['lines.linewidth']    = 2
mpl.rcParams['lines.markersize']   = 6
mpl.rcParams['font.size']          = 18 
mpl.rcParams['legend.fontsize']    = 20
mpl.rcParams['axes.titlesize']     = 18
mpl.rcParams['axes.labelsize']     = 18
mpl.rcParams['xtick.major.size']   = 6
mpl.rcParams['ytick.major.size']   = 6
mpl.rcParams['xtick.major.width']  = 1
mpl.rcParams['ytick.major.width']  = 1
mpl.rcParams['xtick.direction']    = 'in'
mpl.rcParams['ytick.direction']    = 'in'
mpl.rcParams['figure.titlesize']   = 20
mpl.rcParams['figure.figsize']     = [8,3]

#scipy.sparse.linalg.expm

# Set parameter
L  = 10
J = 1
deltaI = 0.5
deltaF = 5.0


hz = 0
n_lowest_eigenvalues = 1


# Step 1: creating initial hamiltonian
rows, cols, data = Hamiltonian.get_hamiltonian_sparse(L, J, deltaI, hz)
hamiltonian = sp.sparse.csc_matrix((data, (rows, cols)))


# Step 2: evaluate initial GS 
energy, groundstate = sp.sparse.linalg.eigsh(hamiltonian, 
                                   k=n_lowest_eigenvalues,
                                   which='SA',
                                   return_eigenvectors=True,
                                   maxiter=1000)

groundstate = groundstate[:,0]
print("ground state energy (per site):", energy/L, flush=True)

plt.plot(abs(groundstate)**2)

# Step 3: creating final hamiltonian
rows, cols, data = Hamiltonian.get_hamiltonian_sparse(L, J, deltaF, hz)
hamiltonian = sp.sparse.csc_matrix((data, (rows, cols)))

start=time.time()
en, vec = sp.linalg.eigh(hamiltonian.toarray())
print('time : ', time.time()-start)
del(hamiltonian)
cn = groundstate@vec
#del(vec)

time = np.linspace(0,5,1000,endpoint=False)
LE = [] 
rate = []
f1 = open('LE_Fig1a_N12.txt','w')
for t in time:
    LE.append(abs(np.sum(np.abs(cn)**2*np.exp(-1j*en*t)))**2)
    rate.append(-np.log(LE[-1])/L)
    f1.write('%f  %f\n'%(t,rate[-1]))
    
f1.close()

plt.figure()
plt.plot(time,rate,'-r')
#plt.plot(time,LE,'-b')
plt.grid('on')
#plt.ylim([0,1])
plt.xlim([0,5])  
    
print("Run successful")





    




