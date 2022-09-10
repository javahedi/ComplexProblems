#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 10:01:41 2019

@author: javad
"""

import random
import numpy as np
import math

def get_hamiltonian_sparse(L, N, J0, hz, alpha,eta,bc='obc'):
    '''
    Creates the Hamiltonian of the long-ranged,
    decaying with a power law XX model
    on a linear chain lattice with periodic boundary conditions.

    Args:
        L(int): length of chain
        N(int): number of spin

        J0 (float): coupling constant for XX
        alpha (float): power exponent of coupling constant

    Returns:
        (hamiltonian_rows, hamiltonian_cols, hamiltonian_data) where:
        hamiltonian_rows (list of ints): row index of non-zero elements
        hamiltonian_cols (list of ints): column index of non-zero elements
        hamiltonian_data (list of floats): value of non-zero elements
    '''

    def get_site_value(state, site):
        ''' Function to get local value at a given site '''
        return (state >> site) & 1
    
    def get_bit_value(state, site):
         ''' Function to get local value at a given site '''
         return ((state&(1<<site))!=0)

    def hilbertspace_dimension(N):
        ''' return dimension of hilbertspace '''
        return 2**N
     
    # Define chain lattice
    #heisenberg_bonds = [(site, (site+1)%N) for site in range(N)]
    #np.random.seed()
    #site_index = np.random.randint(L,size=N)
    #site_index = random.sample(range(L), N)
    #site_index.sort()
    site_index = [13,20,41,42,50,66,73,76,85,97,100,110,115,120] 
    #site_index =  [14,16,25,29,31,48,49,63,66,73,76,82,110,112,119,130]
    
    
    xx_bonds = []
    for k in range(1,N-1):  # k-th neighbour
        Nmax = N if bc == 'pbc' else N-k
        for i in range(Nmax):
             xx_bound.append((i,np.mod(i+k,N)))
\
                
    
    # Empty lists for sparse matrix
    hamiltonian_rows = []
    hamiltonian_cols = []
    hamiltonian_data = []

    # Run through all spin configurations
    for state in range(hilbertspace_dimension(N)):
        
        # Apply Ising bonds
        ising_diagonal = 0
        for bond in xx_bonds:
            if get_site_value(state, bond[0]) == get_site_value(state, bond[1]):
                ising_diagonal += 0/4 
            else:
                ising_diagonal -= 0/4 
        hamiltonian_rows.append(state)
        hamiltonian_cols.append(state)
        hamiltonian_data.append(ising_diagonal)
        
        

        # Apply XX exchange interaction
        for bond in xx_bonds:           
            n,m=bond[0],bond[1]
            z=abs(site_index[n]-site_index[m])
            J = np.exp(-z/eta)/(z**alpha)
            #J = J0/(z**alpha)
            
            
            flipmask = (1 << bond[0]) | (1 << bond[1])
            if get_site_value(state, bond[0]) != get_site_value(state, bond[1]):
                new_state = state ^ flipmask
                hamiltonian_rows.append(new_state)
                hamiltonian_cols.append(state)
                hamiltonian_data.append(J/2)

        

    return hamiltonian_rows, hamiltonian_cols, hamiltonian_data
