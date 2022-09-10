#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 09:57:38 2019

@author: javadvahedi
"""

import numpy as np
import scipy as sp
#np.set_printoptions(precision=4)
#from scipy import sparse
#from scipy.sparse import linalg

#paper arXiv:1601.07458   Eq.31
def get_ptrace(rho_in,da,db):
    rho_out=np.zeros((da,da))
    for j in range(da):
        for k in range(j,da):
            for l in range(db):
                rho_out[j,k] += rho_in[j*db+l,k*db+l]
            if j!=k: 
                rho_out[k,j] += np.conjugate(rho_out[j,k])
    return rho_out
    #rho_out[rho_out<1e-20]=0
    #return sp.sparse.csr_matrix(rho_out)











'''
def ptrace(A,x):
    reshaped = A.reshape(2,2,2,2)
    if x==0:
        return np.einsum('ijik->jk',reshaped)
    else:
        return np.einsum('jiki->jk',reshaped)
    
    
a=np.array([[1,0],[0,0]])
b=np.array([[0,0],[0,1]])
S = np.kron(a,b)
#S = sp.sparse.csr_matrix(np.kron(a,b))

print(S)
print(a)
print(b)
print(ptrace(S,0))
print(ptrace(S,1))



a=np.arange(25).reshape(5,5)
np.einsum('ii',a)==np.trace(a)
np.einsum('ii->i',a)==np.diag(a)
np.einsum('ij->i',a)==np.sum(a,axis=1)
np.einsum('ij->j',a)==np.sum(a,axis=0)
'''

