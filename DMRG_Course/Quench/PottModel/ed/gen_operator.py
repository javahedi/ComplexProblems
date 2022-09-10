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
    
    omega = np.exp(-1j*2*np.pi/p)

    ## local operator ##
    data = np.ones  (p)
    row  = np.arange(p)
    col  = [  2,  0,  1]
    t    = sparse.csr_matrix((data,(row,col)), shape=(p,p))
    
    data = [1.0, omega, omega**2]
    row  = np.arange(p)
    col  = np.arange(p)
    s    = sparse.csr_matrix((data,(row,col)), shape=(p,p), dtype=np.complex128)
    
    data = np.ones(p**L)
    row  = np.arange(p**L)
    col  = np.arange(p**L)
    I    = sparse.csr_matrix((data,(row,col)), shape=(p**L,p**L))
    
    
    identity_list = []
    sigma_list    = []
    tau_list      = [] 
    
    for i_site in range(L):
        if i_site==0: 
            sigma = s 
            tau   = t
        else: 
            sigma = sparse.csr_matrix(np.eye(p)) 
            tau = sparse.csr_matrix(np.eye(p)) 
            
        for j_site in range(1,L): 
            if j_site==i_site: 
                sigma = sparse.kron(sigma, s, 'csr')
                tau   = sparse.kron(tau  , t, 'csr') 
            else: 
                sigma = sparse.kron(sigma , np.eye(p),'csr') 
                tau   = sparse.kron(tau   , np.eye(p),'csr') 
                
        sigma_list.append(sigma)
        tau_list.append(tau) 
        identity_list.append(I)
        
     

    return identity_list,sigma_list,tau_list 

