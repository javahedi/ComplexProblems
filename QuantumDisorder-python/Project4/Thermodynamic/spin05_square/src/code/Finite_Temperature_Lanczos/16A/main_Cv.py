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
L = 16 # system size
m = 50 # dimensio of Krylov space
J1,J2 = 1.0,0.1
hz =0.01
#
N_samples = 10 # of samples to approximate thermal expectation value with
#
#T = np.logspace(-3,3,51,base=10) # temperature vector
T = np.linspace(0.01,2,100) # temperature vector
beta = 1.0/(T+1e-15) # inverse temperature vector
#

#with twist periodic boubdary
J1_bond = [(0,1)  , (1,2)  , (2,3)  , (4,5)  ,
           (5,6)  , (6,7)  , (8,9)  , (9,10) ,
           (10,11), (12,13), (13,14), (14,15),
           (0,4)  , (2,6)  , (5,9)  , (7,11) ,
           (8,12) , (10,14), (0,3)  , (4,7)  ,
           (8,11) , (12,15), (1,13) , (3,15) ]
           
J2_bond = [(0,5)  , (0,2) , (1,4)  , (1,6)  , (1,3)  ,
           (2,5)  , (2,7) , (4,6)  , (4,9)  , (5,7)  ,
           (5,8)  , (5,10), (6,9)  , (6,11) , (7,10) ,
           (8,10) , (8,13), (9,11) , (9,12) , (9,14) ,
           (10,13),(10,15), (11,14), (12,14), (13,15),
           (0,7)  , (0,13), (0,15) , (1,12) , (1,14) ,
           (2,13) , (2,15), (3,12) , (3,14) , (4,3)  ,
           (4,11) , (7,8) , (8,15) , (11,12)         ]
           

bonds = J1_bond + J2_bond


J_list = []
for b in J1_bond:
    J_list.append(0.25 * J1 )
for b in J2_bond:
    J_list.append(0.25 * J2 )



# %% Generate operators
op  =  gen_operator.Operator(L)

hamiltonian = 0
# ZZ
for OP,J in zip(['sx','sx'],[J_list,J_list]):
    op.gen_s0sxsysz(ope=OP)
    ham = gen_hamiltonian.Hamiltonian(op.s_list)
    ham.gen_nn_int(J,bonds)
    hamiltonian +=  ham.H
        
    # create magnetization-squared operator
    #ham.gen_magnetization_squared()
    #M2 = ham.MM
    
    
    

# field in z direction
#op.gen_s0sxsysz(ope='sz')
#ham = gen_hamiltonian.Hamiltonian(op.s_list)
#ham.gen_onsite_field(L*[hz])
#hamiltonian +=  ham.h


n_lowest_eigenvalues = 1
E0 = sp.sparse.linalg.eigsh(hamiltonian,
                                     k=n_lowest_eigenvalues,
                                     which='SA',
                                     return_eigenvectors=False,
                                     maxiter=1E4)


##### finite temperature methods #####
#
# preallocate lists to store results from iterations
#E2_FT_list = []
#E_FT_list = []
#Z_FT_list = []
#Cv_FT_list = []

E2_LT_list = []
E_LT_list = []
Z_LT_list = []
#Cv_LT_list = []
#
# allocate memory for lanczos vectors
out = np.zeros((m,hamiltonian.shape[0]),dtype=np.float64)
#
measurments = {"E":hamiltonian, "E2":hamiltonian**2}
# calculate iterations
for i in range(N_samples):
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
    
    results_LT,Id_LT = LTLM_static_iteration(measurments,E,V,lv,beta=beta)
    
    E_LT_list.append(results_LT["E"])
    E2_LT_list.append(results_LT["E2"])
    Z_LT_list.append(Id_LT)
    
    
    print(f'loop {i}, time : {time.time()-start}')
#
# calculating error bars on the expectation values
E_LT ,dE_LT  = bootstrap_mean(E_LT_list,Z_LT_list)
E2_LT,dE2_LT = bootstrap_mean(E2_LT_list, Z_LT_list)

cv = (E2_LT -E_LT**2)/T**2

np.savetxt('thermodynamic_FTLM.txt',np.column_stack([T,E_LT,E2_LT,cv]))

##### plot results #####
#
# setting up plot and inset
h=4.2 # figure aspect ratio parameter
f,ax = plt.subplots(figsize=(1.5*h,h))
ax.plot(T,cv,'or',label="FTLM")
ax.plot(T,cv,'-k',label="FTLM")
   
plt.show()

