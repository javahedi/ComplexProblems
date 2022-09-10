#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 12:26:14 2020

@author: javadvahedi
"""

import numpy as np
import scipy
import scipy.sparse as sparse



def gen_nn_int(op_list,op_list2, J_list, f_list , phi_list, bc='obc'):
    return gen_interaction_kdist(op_list, op_list2, J_list, f_list, phi_list, 1, bc)


# generates \sum_i O_i O_{i+k} type interactions
def gen_interaction_kdist(op_list, op_list2, J_list, f_list, phi_list, k=1, bc='obc'):
    L= len(op_list)
    H = sparse.csr_matrix(op_list[0].shape)
    Lmax = L if bc == 'pbc' else L-k
    for i in range(Lmax):
        j = np.mod(i+k,L)
        H += J_list[i] * np.exp( 1j*phi_list[i]) * op_list[i]        @ op_list[j].conj()
        H += J_list[i] * np.exp(-1j*phi_list[i]) * op_list[i].conj() @ op_list[j]
        
    for i in range(L):   
        H += f_list[i] * op_list2[i]
        H += f_list[i] * op_list2[i].T 
           
    return H



