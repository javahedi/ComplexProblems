#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 15:32:53 2019

@author: javadvahedi
"""

import numpy as np
import scipy as sp
from scipy import sparse


def get_trace(Psi0,da,db):
    
    '''
    rho = np.zeros((da,da))
    for i in range(da):
        for j in range(i,da):
            for m in range(db):   
                rho[i,j] += Psi0[m+i*db]*Psi0[m+i*db+(j-i)*db]
            if i!=j: 
                rho[j,i] += np.conjugate(rho[i,j])
                
    #rho[rho<1e-10]=0
    return rho
    
    '''

    rows = []
    cols = []
    data = []
    for i in range(da):
        for j in range(da):
            s=0
            for m in range(db):
                #print(i,j,m+i*db,m+i*db+(j-i)*db)
                s += Psi0[m+i*db]*Psi0[m+i*db+(j-i)*db]
                #s += Psi0[i*db+m]*Psi0[j*db+m]
                
                
                
            rows.append(i)
            cols.append(j)
            data.append(s)
            
    
    rho = sp.sparse.csr_matrix((data, (rows, cols)))
    #rho[rho<+1e-10]=0
    return rho
    
            
    
    