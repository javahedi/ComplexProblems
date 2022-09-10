#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on  Dec  20  2020

@author: javadvahedi
"""

import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
from scipy.sparse import linalg
import matplotlib.pyplot as plt
import time
import gen_operator as GO
import gen_op_total as GOT
import gen_nn_int as Gint
import tqdm

import sys,os
lanczos_path = os.path.join(os.getcwd(),"lanczos")
sys.path.insert(0,lanczos_path)
from lanczos import lanczos_full,lanczos_iter,lin_comb_Q_T,expm_lanczos



###########################
# %% system parametrs
###########################

L = 8#int(np.loadtxt('ChainLength.txt'))
p = 3#int(np.loadtxt('Flavour.txt'))
# %% Generate operators
I0, sigma, tau = GO.gen_operator(L,p)
print(" Operators are created ",              flush=True)
print(" Hilber dim : "         , I0[0].shape, flush=True)


#%%
# Step 1: creating initial hamiltonian
J_list   = L * [-1.0]
f_list   = L * [-3.]
phi_list = L * [0.0]
hamiltonian  = Gint.gen_nn_int(sigma,tau,J_list, f_list, phi_list,bc='pbc')


# Step 2: evaluate initial GS 
energy, groundstate = sp.sparse.linalg.eigsh(hamiltonian,
                                             ncv=20,
                                             k=1,
                                             which='SA',
                                             return_eigenvectors=True,
                                             maxiter=1E4)

groundstate = groundstate[:,0]
print("ground state energy (per site):", energy[0]/L, flush=True)
J_list   = L * [-1.0]
f_list   = L * [-0.4]
phi_list = L * [0.0]
hamiltonian  = Gint.gen_nn_int(sigma,tau,J_list, f_list, phi_list,bc='pbc')


################################
#  #excat result
# from PM to FM
J = -1
m = 100
t = np.linspace(0,4,m,endpoint=False)
dt  = t[1]-t[0]

g1 = np.exp(-1j*J*t)*(np.exp(-3*1j*J*t)+2)/3
g2 = np.exp(-1j*J*t)*(np.exp(-3*1j*J*t)-1)/3
Gt = g1**L + 2 * g2**L
Rt = -np.log(np.abs(Gt)**2)/L
################################




LE = []
rt = []
f1 = open('data_L_'+str(L)+'.txt','w')
psi     = groundstate
expH_v0 = groundstate
for i in tqdm.tqdm(range(m)):
    E, V, Q_T = lanczos_iter(hamiltonian,expH_v0,20)
    expH_v0 = expm_lanczos(E,V,Q_T,a=-1j*dt)
    LE.append(abs(expH_v0.dot(groundstate))**2)
    rt.append(-np.log(LE[-1])/L)
    f1.write('%f    %f   %f   %f\n'%(t[i],LE[-1],rt[-1],Rt[i]))
    f1.flush()


f1.close()


print("Run successful")






