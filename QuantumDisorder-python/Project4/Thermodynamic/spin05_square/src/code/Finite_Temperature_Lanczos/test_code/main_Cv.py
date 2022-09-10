#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import print_function, division



import sys,os
lanczos_path = os.path.join(os.getcwd(),"lanczos")
sys.path.insert(0,lanczos_path)
from lanczos import lanczos_full,lanczos_iter,FTLM_static_iteration,LTLM_static_iteration

from mpl_toolkits.axes_grid1.inset_locator import inset_axes


import matplotlib.pyplot as plt
import math
import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
from scipy.sparse import linalg
import gen_operator
import gen_hamiltonian 
import time
import random
import gen_diagonal_ME as GME



#
np.random.seed(1203901) # fix seed
#
#
def bootstrap_mean(O_r,Id_r,n_bootstrap=100):
    """
    Uses boostraping to esimate the error due to sampling.

    O_r: numerator
    Id_r: denominator
    n_bootstrap: bootstrap sample size

    """
    O_r = np.asarray(O_r)
    Id_r = np.asarray(Id_r)
    #
    avg = np.nanmean(O_r,axis=0)/np.nanmean(Id_r,axis=0)
    n_Id = Id_r.shape[0]
    #n_N = O_r.shape[0]
    #
    i_iter = (np.random.randint(n_Id,size=n_Id) for i in range(n_bootstrap))
    #
    bootstrap_iter = (np.nanmean(O_r[i,...],axis=0)/np.nanmean(Id_r[i,...],axis=0) for i in i_iter)
    diff_iter = ((bootstrap-avg)**2 for bootstrap in bootstrap_iter)
    err = np.sqrt(sum(diff_iter)/n_bootstrap)
    #
    return avg,err




#
#
##### define system parameters #####
#
L = 12 # system size
m = 50 # dimensio of Krylov space
J1,J2 = 1.0,1.0
hz =0.01
#
N_samples = 100 # of samples to approximate thermal expectation value with
#
#T = np.logspace(-3,3,51,base=10) # temperature vector
T = np.linspace(0.1,5,10) # temperature vector
beta = 1.0/(T+1e-15) # inverse temperature vector
#

# Define chain lattice

#ladder with periodic boubdary
bonds = [(0,1), (1,2), (2,3), (3,4 ), (4,5),
         (6,7), (7,8), (8,9), (9,10), (10,11),
         (0,6), (1,7), (2,8), (3,9 ), (4,10), (5,11),
         (0,5),(6,11)]

J_list = []
for bond in bonds:
    if abs(bond[0]-bond[1])==5:
        J_list.append(0.25 * J1)
    else:
        J_list.append(0.25 * J2)


Jxx = J_list
Jyy = J_list
Jzz = J_list

# %% Generate operators
op  =  gen_operator.Operator(L)

hamiltonian = 0
for OP,J in zip(['sx','sx','sz'],[Jxx,Jxx,Jzz]):
    op.gen_s0sxsysz(ope=OP)
    ham = gen_hamiltonian.Hamiltonian(op.s_list)
    ham.gen_nn_int(J,bonds)
    hamiltonian +=  ham.H
        
    # create magnetization operator
    if OP=='sz':
        ham.gen_onsite_field(L*[0.5/L])
        Mz  = ham.h
        Mz2 = Mz @ Mz

    
    
    

# field in z direction
op.gen_s0sxsysz(ope='sz')
ham = gen_hamiltonian.Hamiltonian(op.s_list)
ham.gen_onsite_field(L*[hz])
hamiltonian +=  ham.h



n_lowest_eigenvalues = 1
E0 = sp.sparse.linalg.eigsh(hamiltonian,
                                     k=n_lowest_eigenvalues,
                                     which='SA',
                                     return_eigenvectors=False,
                                     maxiter=1E4)


##### finite temperature methods #####
#
# preallocate lists to store results from iterations


E2_LT_list = []
E_LT_list  = []
M2_LT_list = []
M_LT_list  = []
Z_LT_list  = []
#
# allocate memory for lanczos vectors
out = np.zeros((m,hamiltonian.shape[0]),dtype=np.float64)
#

measurement ={"E" :hamiltonian,    "M" :Mz,
              "E2":hamiltonian**2, "M2":Mz2 }

# calculate iterations
for i in range(N_samples):
    print(i)
    # generate normalized random vector
    r = np.random.normal(0,1,size=hamiltonian.shape[0])
    r /= np.linalg.norm(r)
    # get lanczos basis
    E,V,lv = lanczos_full(hamiltonian,r,m,eps=1e-8,full_ortho=True)
    # E,V,lv = lanczos_full(hamiltonian,r,m,eps=1e-8,full_ortho=False)
    # E,V,lv = lanczos_iter(hamiltonian,r,m,eps=1e-8)
    # shift energy to avoid overflows
    E -= E0
    # calculate iteration
    #results_FT,Id_FT = FTLM_static_iteration({"E2":hamiltonian**2},E,V,lv,beta=beta)
    results_LT,Id_LT = LTLM_static_iteration(measurement,E,V,lv,beta=beta)
    # save results to a list
    E_LT_list.append(results_LT["E"])
    M_LT_list.append(results_LT["M"])
    E2_LT_list.append(results_LT["E2"])
    M2_LT_list.append(results_LT["M2"])
    Z_LT_list.append(Id_LT)
    
    
#
# calculating error bars on the expectation values
E_LT ,dE_LT  = bootstrap_mean(E_LT_list,  Z_LT_list)
E2_LT,dE2_LT = bootstrap_mean(E2_LT_list, Z_LT_list)

M_LT ,dM_LT  = bootstrap_mean(M_LT_list,  Z_LT_list)
M2_LT,dM2_LT = bootstrap_mean(M2_LT_list, Z_LT_list)

cv = (E2_LT -E_LT**2)/T**2
chi = (M2_LT -M_LT**2)/T



##### calculating exact results from full diagonalization #####
#
M2_ED_list = []
M_ED_list  = []

dim_cutoff=5000 # Hilbert space dimension cutoff
if hamiltonian.shape[0] < dim_cutoff: # Hilbert space is not too big to diagonalize on a laptop
    
    #
    # full diagonaization of H
    E,V = np.linalg.eigh(hamiltonian.todense())
    # shift energy to avoid overflows
    E -= E[0]
    
    for b in beta:
        print(b)
        # get boltzmann weights for each temperature
        W = np.exp(-E*b)
        
        O  = GME.gen_diagonal_ME(Mz , V)
        O2 = GME.gen_diagonal_ME(Mz2, V)
        # calculate trace
        M_ED_list.append(np.sum(O@W)/np.sum(W))
        M2_ED_list.append(np.sum(O2@W)/np.sum(W))


M_ED = np.asarray(M_ED_list)
M2_ED = np.asarray(M_ED_list)
#
#

Temp, Z,  E, E2, Cv = np.loadtxt('thermodynamic_FTLM.txt',skiprows=2, unpack=True)
##### plot results #####
#
# setting up plot and inset
h=4.2 # figure aspect ratio parameter
f ,  ax = plt.subplots(figsize=(1.5*h,h))
axinset = inset_axes(ax, width="45%", height="45%", loc="upper right")
#axs = [ax,axinset]
#
# plot results for FTLM and LTLM.

#axs[0].errorbar(T,E_LT ,dE_LT , marker=".",label="<E>",zorder=-1)
#axs[0].errorbar(T,E2_LT,dE2_LT, marker=".",label="<E2>", zorder=-2)
    
#axs[0].plot(T,E_LT,label="FTLM")
#axs[0].plot(T,E,label="ED")

#axs[1].loglog(T,abs(E_LT-E))
ax.plot(T,cv,'or',label="FTLM")
ax.plot(Temp,Cv,'--b',label="ED")

axinset.plot(T,M_LT**2,'or')
axinset.plot(T,M2_LT,'sg')
axinset.plot(T,M_ED**2,'*r')
axinset.plot(T,M2_ED,'^g')
#axinset.plot(T,chi,'.b')
   
plt.legend()
plt.show()

