#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 10:01:41 2019

@author: javad
"""

import random
#import numpy as np
import math

def get_hamiltonian_sparse(L, J, delta, hz):
    '''
    Creates the Hamiltonian of the long-ranged,
    decaying with a power law XX model
    on a linear chain lattice with periodic boundary conditions.

    Args:
        L(int): length of chain
        

        J (float): coupling constant for XX+YY
        delta (float): coupling constant for ZZ

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

    def hilbertspace_dimension(L):
        ''' return dimension of hilbertspace '''
        return 2**L
     
    # Define chain lattice
    heisenberg_bonds = [(site, (site+1)%L) for site in range(L)]

   
  
    
    # Empty lists for sparse matrix
    hamiltonian_rows = []
    hamiltonian_cols = []
    hamiltonian_data = []

    # Run through all spin configurations
    for state in range(hilbertspace_dimension(L)):
        
        # Apply ZZ bonds
        ising_diagonal = 0
        for bond in heisenberg_bonds:
            if get_site_value(state, bond[0]) == get_site_value(state, bond[1]):
                ising_diagonal += J/4 + hz/2
            else:
                ising_diagonal -= J/4 - hz/2
        hamiltonian_rows.append(state)
        hamiltonian_cols.append(state)
        hamiltonian_data.append(ising_diagonal)
        
        
        
        # Apply XX exchange interaction
        for bond in heisenberg_bonds:           
            flipmask = (1 << bond[0]) | (1 << bond[1])
            if get_site_value(state, bond[0]) != get_site_value(state, bond[1]):
                new_state = state ^ flipmask
                hamiltonian_rows.append(new_state)
                hamiltonian_cols.append(state)
                hamiltonian_data.append((1+delta)*J/4)
                
                
        '''
        # Apply transverse field
        for site in range(L):
            # Flip spin at site
            new_state = state ^ (1 << site)
            hamiltonian_rows.append(new_state)
            hamiltonian_cols.append(state)
            hamiltonian_data.append(hx/2)
        '''

        

    return hamiltonian_rows, hamiltonian_cols, hamiltonian_data
