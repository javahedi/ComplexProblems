#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March-06  2021
@author: javadvahedi
"""
import math
import sys,os
#os.environ['OMP_NUM_THREADS']='1' # set number of OpenMP threads to run in parallel
os.environ['MKL_NUM_THREADS']='2' # set number of MKL threads to run in parallel

import matplotlib.pylab  as plt
import numpy             as np
from scipy.linalg        import eigh
from scipy.sparse.linalg import eigsh
from gen_diagprojector   import gen_diagprojector
from gen_op_total        import gen_op_total
from gen_diagonal_ME     import gen_diagonal_ME
from Entanglement        import SvN, Mutual
import gen_operator
import gen_hamiltonian
import Time_Mem
import time
import random




def run(N,K):

    #--------------------------------
    # %%  system parametrs
    J                    = 1.0
    projection           = False
    full_diagonalization = True
    groun_state          = False
    bc                   ='pbc'

    #--------------------------------
    # %% Define  lattice bounds and couplings

    Nmax = N if bc == 'pbc' else N-1
    bonds = [(i,np.mod(i+1,N)) for i in range(Nmax)]
    Jxx = [0.5 * J] * len(bonds)
    Jyy = Jxx
    Jzz = [J] * len(bonds)
   

    #--------------------------------
    # %% Generate operators
    op  =  gen_operator.Operator(N)

    hamiltonian = 0.

    #--------------------------------
    #%% XXZ term
    for OP,J in zip(['sx','sy','sz'],[Jxx,Jyy,Jzz]):
        op.gen_s0sxsysz(ope=OP)
        ham = gen_hamiltonian.Hamiltonian(op.s_list)
        ham.gen_nn_int(J,bonds)
        hamiltonian +=  ham.H
        
        
    

    #--------------------------------
    #%% field term
    #hx = [0.001] * N
    #op.gen_s0sxsysz(ope='sx')
    #ham = gen_hamiltonian.Hamiltonian(op.s_list)
    #ham.gen_onsite_field(hx)
    #hamiltonian +=  ham.field
    
    #%% single-site anisotropy term
    Kz = [-K] * N
    op.gen_s0sxsysz(ope='sz')
    ham = gen_hamiltonian.Hamiltonian(op.s_list)
    ham.gen_onsite_anisotropy(Kz)
    hamiltonian +=  ham.aniotropy


    #print( ' Before projection : ', flush=True)
    #print( ' Hamiltonian   dim : ', hamiltonian.shape, flush=True)
    #print('%'*40)




    #--------------------------------
    # %% projection
    if projection:
        val    = 0
        op.gen_s0sxsysz(ope='sz')
        Sztot  = gen_op_total(op.s_list)
        
        P     = gen_diagprojector(Sztot.diagonal(), val)
        hamiltonian  = P @ hamiltonian @ P.T
        #print( ' After projection : ', flush=True)
        #print( ' Hamiltonian   dim : ', hamiltonian.shape, flush=True)
        #print('%'*40)


    #--------------------------------
    # %% finding energy
    if full_diagonalization:
        eigs = eigh(hamiltonian.todense(), eigvals_only=True)
        np.savetxt(f'energy_L{N}.txt',eigs)
        
    if groun_state:
        eigs = eigsh(hamiltonian, k=2,
                                  which='SA',
                                  return_eigenvectors=False,
                                  maxiter=1E4)
                                  
    return eigs



if __name__ == "__main__":
    #--------------------------------
    start = time.time()
    N=10
    K=0.0
    eigs = run(N,K)
    #Temp = np.linspace(0.1,5,100)
    Temp = np.logspace(-2,1,100,base=10)
    Cv = []
    f1 = open(f'thermodynamic_ed_S1_L{N}_K_{K}.txt','w')
    f1.write('Temperature      Z      <E>     <E^2>    Cv\n')
    f1.write('============================================\n')
    
    for t in Temp:
       Z = np.sum(np.exp(-eigs/t))
       ave_E  = np.sum(eigs    *np.exp(-eigs/t))/ Z
       ave_E2 = np.sum(eigs**2 *np.exp(-eigs/t))/Z
       Cv.append((ave_E2-ave_E**2)/t**2)
       f1.write('%f   %f    %f    %f    %f\n'%(t,Z, ave_E, ave_E2, Cv[-1]))
       f1.flush()
    
    f1.close()
    plt.plot(Temp,Cv,'.-b')
    plt.show()
    #--------------------------------
    #%%  memory check
    _=Time_Mem.check(start,date=True)
    
    
