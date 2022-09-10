#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 12:21:23 2020

@author: javadvahedi
"""
import scipy
import scipy.sparse as sparse
import scipy.sparse.linalg as spalin
import numpy as np

def gen_operator(L,p):
    
    ## local operator ##
    data = np.ones(p-1)
    row = np.arange(p-1)
    col = np.arange(1,p)
    f   = sparse.csr_matrix((data,(row,col)),shape=(p,p))
    
    data = np.arange(p)
    row  = np.arange(p)
    col  = np.arange(p)
    n    = sparse.csr_matrix((data,(row,col)),shape=(p,p))
    
    data = np.ones(p**L)
    row  = np.arange(p**L)
    col  = np.arange(p**L)
    I = sparse.csr_matrix((data,(row,col)),shape=(p**L,p**L))
    
    
    f0_list = []
    fu_list = []
    n_list  = [] 
    
    for i_site in range(L):
        if i_site==0: 
            F = f 
            N = n
        else: 
            F = sparse.csr_matrix(np.eye(p)) 
            N = sparse.csr_matrix(np.eye(p)) 
            
        for j_site in range(1,L): 
            if j_site==i_site: 
                F = sparse.kron(F,f, 'csr')
                N = sparse.kron(N,n, 'csr') 
            else: 
                F = sparse.kron(F,np.eye(p),'csr') 
                N = sparse.kron(N,np.eye(p),'csr') 
                
        fu_list.append(F)
        n_list.append(N) 
        f0_list.append(I)

    return f0_list,fu_list,n_list 
