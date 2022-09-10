# -*- coding: utf-8 -*-
"""
Function to create the Hamiltonian of the spin 1/2 Heisenberg model
with staggered magnetic field, Sz conservation implemented

:author: Alexander Wietek
:email: alexander.wietek@uibk.ac.at
:year: 2018
"""
from __future__ import absolute_import, division, print_function
import random
def get_hamiltonian_sparse(L, N, J0, hz, alpha, sz):
    '''
    Creates the Hamiltonian of the Heisenberg model in a staggered magnetic
    field on a linear chain lattice with periodic boundary conditions.

    Args:
        J0 (float): coupling constant for XX
        alpha (float): power exponent of coupling constant
        hz (float): coupling constant for field
        sz (int): total Sz

    Returns:
        (hamiltonian_rows, hamiltonian_cols, hamiltonian_data) where:
        hamiltonian_rows (list of ints): row index of non-zero elements
        hamiltonian_cols (list of ints): column index of non-zero elements
        hamiltonian_data (list of floats): value of non-zero elements
    '''

    def get_site_value(state, site):
        ''' Function to get local value at a given site '''
        return (state >> site) & 1

    # Functions to create sz basis
    def first_state(N, sz):
        ''' Return first state of Hilbert space in lexicographic order '''
        n_upspins = N//2 + sz
        return (1 << n_upspins) - 1

    def next_state(state):
        '''
        Return next state of Hilbert space in lexicographic order

        This function implements is a nice trick for spin 1/2 only,
        see http://graphics.stanford.edu/~seander/bithacks.html
        # Next Bit Permutation for details
        '''
        t = (state | (state - 1)) + 1
        return t | ((((t & -t) // (state & -state)) >> 1) - 1)

    def last_state(L, sz):
        ''' Return last state of Hilbert space in lexicographic order '''
        n_upspins = N//2 + sz
        return ((1 << n_upspins) - 1) << (N - n_upspins)


    # check if sz is valid
    assert (sz <= (N // 2 + N % 2)) and (sz >= -N//2)

    # Create list of states with fixed sz
    basis_states = []
    state = first_state(N, sz)
    end_state = last_state(N, sz)
    while state <= end_state:
        basis_states.append(state)
        state = next_state(state)




    # Define chain lattice
    #heisenberg_bonds = [(site, (site+1)%N) for site in range(N)]
    site_index = random.sample(range(L), N)
    site_index.sort()
    
    
    
    xx_bonds = []
    for i in range(N):
        for j in range(i+1,N):  # N-1  for PBC
                xx_bonds.append((i,j))
    



    # Empty lists for sparse matrix
    hamiltonian_rows = []
    hamiltonian_cols = []
    hamiltonian_data = []

    # Run through all spin configurations
    for state_index, state in enumerate(basis_states):

        # Apply diagonal Ising bonds
        diagonal = 0
        for bond in xx_bonds:
            if get_site_value(state, bond[0]) == get_site_value(state, bond[1]):
                diagonal += 0/4 + hz/2
            else:
                diagonal -= 0/4 + hz/2



        hamiltonian_rows.append(state_index)
        hamiltonian_cols.append(state_index)
        hamiltonian_data.append(diagonal)

        # Apply exchange interaction
        for bond in xx_bonds:
            n,m=bond[0],bond[1]
            J = J0 * abs(site_index[n]-site_index[m])**-alpha
            flipmask = (1 << bond[0]) | (1 << bond[1])
            if get_site_value(state, bond[0]) != get_site_value(state, bond[1]):
                new_state = state ^ flipmask
                new_state_index = basis_states.index(new_state)
                hamiltonian_rows.append(state_index)
                hamiltonian_cols.append(new_state_index)
                hamiltonian_data.append(J/2)

    return hamiltonian_rows, hamiltonian_cols, hamiltonian_data, site_index
