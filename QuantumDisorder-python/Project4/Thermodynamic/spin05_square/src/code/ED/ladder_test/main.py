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




def run():

    #--------------------------------
    # %%  system parametrs
    N                    = 12
    J1,J2                = 1.0,1.0
    projection           = False
    bc                   ='obc'
    full_diagonalization = True
    groun_state          = False


    #--------------------------------
    # %% Define  lattice bounds and couplings
    
    '''
    bonds = []
    for k in range(1,3):
        Nmax = N if bc == 'pbc' else N-k
        for i in range(Nmax):
             bonds.append((i,np.mod(i+k,N)))
             #print(k,bonds[-1])

    J_list = []
    for bond in bonds:
        if abs(bond[0]-bond[1])==1:
            J_list.append(0.25 * J1)
        if abs(bond[0]-bond[1])==2:
            J_list.append(0.25 * J2)
        
    '''
    
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
    
    
    Jxx = J_list
    Jyy = J_list
    Jzz = J_list



    #--------------------------------
    # %% Generate operators
    op  =  gen_operator.Operator(N)

    hamiltonian = 0

    #--------------------------------
    #%% XXZ term
    for OP,J in zip(['sx','sx','sz'],[Jxx,Jxx,Jzz]):
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
    
    eigs = run()
    Temp = np.linspace(0.1,5,100)
    Cv = []
    f1 = open('thermodynamic.txt','w')
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
    
    
