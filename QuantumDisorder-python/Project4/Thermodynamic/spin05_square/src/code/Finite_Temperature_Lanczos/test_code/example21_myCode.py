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
J=-0.5
hx=0.5
#
N_samples = 100 # of samples to approximate thermal expectation value with
#
T = np.logspace(-3,3,51,base=10) # temperature vector
beta = 1.0/(T+1e-15) # inverse temperature vector

#

# Define chain lattice
bonds = [(i,(i+1)%L) for i in range(L)]
J_list = [] 
for bond in bonds:
    J_list.append(J)

"""   
#ladder with periodic boubdary
bonds = [(0,1), (1,2), (2,3), (3,4), (4,5),
         (6,7), (7,8), (8,9), (9,10), (10,11),
         (0,6), (1,7), (2,8), (3,9), (4,10), (5,11),
         (0,5),(6,11)]

J_list = []
for bond in bonds:
    if abs(bond[0]-bond[1])==5:
        J_list.append(0.25 * J1)
    else:
        J_list.append(0.25 * J2)
"""

Jxx = J_list
Jyy = J_list
Jzz = J_list

# %% Generate operators
op  =  gen_operator.Operator(L)

hamiltonian = 0
# ZZ
for OP,J in zip(['sz'],[Jzz]):
    op.gen_s0sxsysz(ope=OP)
    ham = gen_hamiltonian.Hamiltonian(op.s_list)
    ham.gen_nn_int(J,bonds)
    hamiltonian +=  ham.H
        
    # create magnetization-squared operator
    #ham.gen_magnetization_squared()
    #M2 = ham.MM
    
    
# field in x direction
op.gen_s0sxsysz(ope='sx')
ham = gen_hamiltonian.Hamiltonian(op.s_list)
ham.gen_onsite_field(L*[hx])
hamiltonian +=  ham.h


# field in z direction
op.gen_s0sxsysz(ope='sz')
ham = gen_hamiltonian.Hamiltonian(op.s_list)
ham.gen_onsite_field(L*[1.0/L])


# find ground state
n_lowest_eigenvalues = 1
E0 = sp.sparse.linalg.eigsh(hamiltonian,
                                     k=n_lowest_eigenvalues,
                                     which='SA',
                                     return_eigenvectors=False,
                                     maxiter=1E4)


##### finite temperature methods #####
#
# preallocate lists to store results from iterations
M2_FT_list = []
M2_LT_list = []
Z_FT_list  = []
Z_LT_list  = []
#
# allocate memory for lanczos vectors
out = np.zeros((m,hamiltonian.shape[0]),dtype=np.float64)
#

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
    results_FT,Id_FT = FTLM_static_iteration({"M2":ham.h**2},E,V,lv,beta=beta)
    results_LT,Id_LT = LTLM_static_iteration({"M2":ham.h**2},E,V,lv,beta=beta)

    # save results to a list
    M2_FT_list.append(results_FT["M2"])
    Z_FT_list.append(Id_FT)
    M2_LT_list.append(results_LT["M2"])
    Z_LT_list.append(Id_LT)
 
    
#
# calculating error bars on the expectation values
m2_FT,dm2_FT = bootstrap_mean(M2_FT_list,Z_FT_list)
m2_LT,dm2_LT = bootstrap_mean(M2_LT_list,Z_LT_list)
   

##### calculating exact results from full diagonalization #####
#
dim_cutoff=2000 # Hilbert space dimension cutoff
if hamiltonian.shape[0] < dim_cutoff: # Hilbert space is not too big to diagonalize on a laptop
    #
    # adding more points for smooth line
    #T_new = np.logspace(np.log10(T.min()),np.log10(T.max()),10*len(T))
    T_new = np.linspace(0.01,5,100) # temperature vector
    beta_new = 1.0/(T_new+1e-15)
    #
    # full diagonaization of H
    E,V = sp.linalg.eigh(hamiltonian.todense())#, eigvals_only=True)
    # shift energy to avoid overflows
    E -= E[0]
    # get boltzmann weights for each temperature
    W = np.exp(-np.outer(E,beta_new))
    # get diagonal matrix elements for trace
    e2 = np.einsum("j...,j->...",W,E**2)/np.einsum("j...->...",W)
    e = np.einsum("j...,j->...",W,E)/np.einsum("j...->...",W)
    
    
    CV = (e2-e**2)/T_new**2
    cv = np.diff(e)
    
    
#

#
#
##### plot results #####
##### plot results #####
#
# setting up plot and inset
h=4.2 # figure aspect ratio parameter
f,ax = plt.subplots(figsize=(1.5*h,h))
axinset = inset_axes(ax, width="45%", height="45%", loc="upper right")
axs = [ax,axinset]
#
# plot results for FTLM and LTLM.
for a in axs:
	a.errorbar(T,m2_LT,dm2_LT,marker="s",label="LTLM",zorder=-1)
	a.errorbar(T,m2_FT,dm2_FT,marker=".",label="FTLM",zorder=-2)
	#
	if hamiltonian.shape[0] < dim_cutoff: # hilbert space is not too big to diagonalize on a laptop
		a.plot(T_new,O,label="exact",zorder=0)
	#
	a.set_xscale("log")
    
plt.legend()
plt.show()
