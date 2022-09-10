#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 12:26:14 2020
@author: javadvahedi
"""

import numpy as np
import scipy
import scipy.sparse as sparse



def gen_nn_int(op_list, op_list2, g_list, m_list ,bc='obc'):
    return gen_interaction_kdist(op_list, op_list2, g_list, m_list, 1, bc)


# generates \sum_i O_i O_{i+k} type interactions
def gen_interaction_kdist(op_list, op_list2, g_list, m_list,k=1, bc='obc'):
    L= len(op_list)
    H = sparse.csr_matrix(op_list[0].shape)
    Lmax = L if bc == 'pbc' else L-k
    for i in range(Lmax):
        j = np.mod(i+k,L)
        H += -(1 - g_list[i]) * op_list[i]   @ op_list[j].T  
        H += -(1 - g_list[i]) * op_list[i].T @ op_list[j] 
        H += -g_list[i] * op_list[i]   @ op_list[i]   @ op_list[j].T  @ op_list[j].T 
        H += -g_list[i] * op_list[i].T @ op_list[i].T @ op_list[j]    @ op_list[j] 
        
        
    if m_list !=[]:
       for i in range(L):
            #H +=  m_list[i] * op_list[i] @ op_list[i].T 
            #H +=  m_list[i] * op_list[i] @ op_list[i] @ op_list[i].T @ op_list[i].T 
            H += m_list[i] * op_list2[i]  
    return H

