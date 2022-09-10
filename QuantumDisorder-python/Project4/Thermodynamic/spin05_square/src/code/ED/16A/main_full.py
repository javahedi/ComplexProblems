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




def run(L, val=0, projection=False):

    #--------------------------------
    # %%  system parametrs
    J1,J2                = 1.0,1.0
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
    '''
    
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

    #--------------------------------
    # %% Generate operators
    op  =  gen_operator.Operator(L)

    hamiltonian = 0

    #--------------------------------
    #%% XXZ term
    for OP,J in zip(['sx','sx'],[J_list,J_list]):
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
        op.gen_s0sxsysz(ope='sz')
        Sztot  = gen_op_total(op.s_list)
        
        P     = gen_diagprojector(Sztot.diagonal(), val)
        hamiltonian  = P @ hamiltonian @ P.T
        #print( ' After projection : ', flush=True)
        print( ' Hamiltonian   dim : ', hamiltonian.shape, flush=True)
        #print('%'*40)


    #--------------------------------
    # %% finding energy
    if full_diagonalization:
        eigs = eigh(hamiltonian.todense(), eigvals_only=True)
        
    if groun_state:
        eigs = eigsh(hamiltonian, k=2,
                                  which='SA',
                                  return_eigenvectors=False,
                                  maxiter=1E4)
                                  
    return eigs






#############################
if __name__ == "__main__":
    #--------------------------------
    start = time.time()
    L = 12
    eigs = run(L, val=0, projection=False)
    np.savetxt(f"energy_L{L}_full.txt",np.sort(eigs))
    
    #
    Temp = np.linspace(0.1,5,100)
    Cv = []
    f1 = open('thermodynamic_full.txt','w')
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
    
    
