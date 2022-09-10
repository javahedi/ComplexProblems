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
L = 10 # system size
m = 50 # dimensio of Krylov space
J = 1.0
K = 1.0

#
N_samples = 100 # of samples to approximate thermal expectation value with
#
T = np.logspace(-2,1,100,base=10) # temperature vector
#T = np.linspace(0.01,2,100) # temperature vector
beta = 1.0/(T+1e-15) # inverse temperature vector
#

bc = 'pbc'
Lmax = L if bc == 'pbc' else L-1
bonds = [(i,np.mod(i+1,L)) for i in range(Lmax)]

Jxx = [0.5 * J] * len(bonds)
Jyy = Jxx
Jzz = [J] * len(bonds)


#--------------------------------
# %% Generate operators
op  =  gen_operator.Operator(L)
hamiltonian = 0.0
#--------------------------------
#%% XXZ term
for OP,J in zip(['sx','sy','sz'],[Jxx,Jyy,Jzz]):
    op.gen_s0sxsysz(ope=OP)
    ham = gen_hamiltonian.Hamiltonian(op.s_list)
    ham.gen_nn_int(J,bonds)
    hamiltonian +=  ham.H

#--------------------------------
#%% field term
hz = [0.01] * L
op.gen_s0sxsysz(ope='sz')
ham = gen_hamiltonian.Hamiltonian(op.s_list)
ham.gen_onsite_field(hz)
hamiltonian +=  ham.field
#--------------------------------
#%% single-site anisotropy term
Kz = [-K] * L
op.gen_s0sxsysz(ope='sz')
ham = gen_hamiltonian.Hamiltonian(op.s_list)
ham.gen_onsite_anisotropy(Kz)
hamiltonian +=  ham.aniotropy



#--------------------------------
#  mesurment operators
#--------------------------------
# the squared magnetization M=(\sum_i S_i^z)/L
# \chi(T) = (<M**2> - <M>**2 )/T
hz = [1.0] * L
op.gen_s0sxsysz(ope='sz')
mz = gen_hamiltonian.Hamiltonian(op.s_list)
mz.gen_onsite_field(hz)
M =  mz.field
del mz
#--------------------------------
# SUCEPTILBILITY
# \chi(T) =  (\sum_i  <S_0^zS_i^z> )/T
s0si = gen_hamiltonian.Hamiltonian(op.s_list)
#--------------------------------



#--------------------------------
n_lowest_eigenvalues = 1
E0 = sp.sparse.linalg.eigsh(hamiltonian,
                                     k=n_lowest_eigenvalues, which='SA',
                                     return_eigenvectors=False,
                                     maxiter=1E4)

#--------------------------------
##### finite temperature methods #####
#
# preallocate lists to store results from iterations
E2_FT_list  = []
E_FT_list   = []
M2_FT_list  = []
M_FT_list   = []
Z_FT_list   = []

S0Si_FT_list = []

E2_LT_list = []
E_LT_list  = []
Z_LT_list  = []

#
# allocate memory for lanczos vectors
out = np.zeros((m,hamiltonian.shape[0]),dtype=np.float64)
#
measurments = {"E":hamiltonian, "E2":hamiltonian**2, "M":M, "M2":M**2}
# calculate iterations
for j in range(N_samples):
    start = time.time()
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
    results_FT,Id_FT = FTLM_static_iteration(measurments,E,V,lv,beta=beta)
    #results_LT,Id_LT = LTLM_static_iteration(measurments,E,V,lv,beta=beta)

    #E_LT_list.append(results_LT["E"])
    #E2_LT_list.append(results_LT["E2"])
    #Z_LT_list.append(Id_LT)

    E_FT_list.append(results_FT["E"])
    E2_FT_list.append(results_FT["E2"])
    M_FT_list.append(results_FT["M"])
    M2_FT_list.append(results_FT["M2"])
    Z_FT_list.append(Id_FT)

    S0Si = 0.0
    for i in range(L):
        s0si.gen_nn_int([1.0],[(0,i)])
        results_FT,Id_FT = FTLM_static_iteration({"S0Si":s0si.H},E,V,lv,beta=beta)
        S0Si +=results_FT["S0Si"]
    S0Si_FT_list.append(S0Si)

    print(f'loop {j}, time : {time.time()-start}',flush=True)
#--------------------------------
#
# calculating error bars on the expectation values
#E_LT ,dE_LT  = bootstrap_mean(E_LT_list,Z_LT_list)
#E2_LT,dE2_LT = bootstrap_mean(E2_LT_list, Z_LT_list)

E_FT ,dE_FT  = bootstrap_mean(E_FT_list,Z_FT_list)
E2_FT,dE2_FT = bootstrap_mean(E2_FT_list, Z_FT_list)
M_FT ,dM_FT  = bootstrap_mean(M_FT_list,Z_FT_list)
M2_FT,dM2_FT = bootstrap_mean(M2_FT_list, Z_FT_list)
chi_FT,dchi_FT = bootstrap_mean(S0Si_FT_list, Z_FT_list)

#cv_LT = (E2_LT -E_LT**2)/T**2
#Z_LT = np.asarray(Z_LT_list)
#Z_LT = np.nanmean(Z_LT,axis=0)
#freeEnergy_LT = -T * np.log(Z_LT)
#entropy_LT =  (E_LT - freeEnergy_LT)/T


cv_FT = (E2_FT -E_FT**2)/T**2
Z_FT = np.asarray(Z_FT_list)
Z_FT = np.nanmean(Z_FT,axis=0)
freeEnergy_FT = -T * np.log(Z_FT)
entropy_FT =  (E_FT - freeEnergy_FT)/T
chi_FT = (M2_FT -M_FT**2)/T
chi2_FT = chi_FT/T
#--------------------------------


#np.savetxt(f'thermodynamic_FTLM_L{L}.txt',np.column_stack([T,freeEnergy, E_LT,E2_LT,cv,entropy]),delimiter=' ',fmt='%f')

f1 = open(f'thermodynamic_FTLM_L{L}_K{K}.txt','w')
f1.write('Temperature,  F    <E>     <E^2>    Cv    entropy   chi    chi2   dE_LT      dE2_LT\n')
f1.write('============================================\n')

#f2 = open(f'thermodynamic_LTLM_L{L}.txt','w')
#f2.write('Temperature,  F    <E>     <E^2>    Cv    entropy    dE_LT      dE2_LT\n')
#f2.write('============================================\n')

#--------------------------------

cte=0
for t in T:
    f1.write('%f   %f   %f    %f    %f   %f   %f    %f   %f   %f\n'%(t, freeEnergy_FT[cte], E_FT[cte], E2_FT[cte], cv_FT[cte], entropy_FT[cte], chi_FT[cte],  chi2_FT[cte], dE_FT[cte], dE2_FT[cte]))
    #f2.write('%f   %f   %f    %f    %f   %f   %f    %f\n'%(t, freeEnergy_LT[cte], E_LT[cte], E2_LT[cte], cv_LT[cte], entropy_LT[cte], dE_LT[cte], dE2_LT[cte]))
    cte+=1
f1.flush()
f1.close()
#--------------------------------
#f2.flush()
#f2.close()
##### plot results #####
#
# setting up plot and inset
#h=4.2 # figure aspect ratio parameter
#f,ax = plt.subplots(figsize=(1.5*h,h))
#ax.plot(T,cv,'or',label="FTLM")
#ax.plot(T,cv,'-k',label="FTLM")
#plt.show()
