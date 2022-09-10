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

import gen_s0sxsysz
import gen_hamiltonian 

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
deltaF = 5
hz = L*[0]



s0,x,y,z = gen_s0sxsysz.gen(L)

xx = gen_hamiltonian.Hamiltonian(x)
yy = gen_hamiltonian.Hamiltonian(y)
zz = gen_hamiltonian.Hamiltonian(z)



# Step 1: creating initial hamiltonian

xx.gen_nn_int(k=1,bc='pbc')
hamiltonian = 0.25* deltaI * xx.H

yy.gen_nn_int(k=1,bc='pbc')
hamiltonian = hamiltonian + 0.25*yy.H

zz.gen_nn_int(k=1,bc='pbc')
hamiltonian = hamiltonian +   0.25*zz.H

zz.gen_onsite_field(hz)
hamiltonian = hamiltonian + 0.5*zz.h
  


# Step 2: evaluate initial GS 
n_lowest_eigenvalues = 1
energy, groundstate = sp.sparse.linalg.eigsh(hamiltonian, 
                                   k=n_lowest_eigenvalues,
                                   which='SA',
                                   return_eigenvectors=True,
                                   maxiter=1000)

groundstate = groundstate[:,0]
print("ground state energy (per site):", energy[0]/L, flush=True)


plt.plot(abs(groundstate)**2)

# Step 3: creating final hamiltonian
xx.gen_nn_int(k=1,bc='pbc')
hamiltonian = 0.25* deltaF * xx.H

yy.gen_nn_int(k=1,bc='pbc')
hamiltonian = hamiltonian + 0.25*yy.H

zz.gen_nn_int(k=1,bc='pbc')
hamiltonian = hamiltonian +   0.25*zz.H

zz.gen_onsite_field(hz)
hamiltonian = hamiltonian + 0.5*zz.h
  


start=time.time()
en, vec = sp.linalg.eigh(hamiltonian.toarray())
print('time : ', time.time()-start)
del(hamiltonian)
cn = groundstate@vec
del(vec)

time = np.linspace(0,10,1000,endpoint=False)
LE = [] 
rate = []
f1 = open('LE_Fig1a_N10.txt','w')
for t in time:
    LE.append(abs(np.sum(np.abs(cn)**2*np.exp(-1j*en*t)))**2)
    rate.append(-np.log(LE[-1])/L)
    f1.write('%f  %f  %f\n'%(t,rate[-1],rate[-1]/4))
    
f1.close()

plt.figure()
plt.plot(time,rate,'-r')
plt.grid('on')
plt.xlim([0,time[-1]])
    
print("Run successful")





    




