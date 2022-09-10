#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 12:37:10 2019

@author: javadvahedi
"""

import numpy as np
from numpy import sqrt

def concurrence(rho):
    
    
    sigmax = np.array([[0,   1], [1,  0]])
    sigmay = np.array([[0, -1j], [1j, 0]])
    sigmaz = np.array([[1,   0], [0, -1]])
    
    
    rho_tilde = np.kron(sigmay,sigmay) @ rho @np.kron(sigmay,sigmay)
    R = rho @ rho_tilde
    eigs = np.linalg.eigvals(R)
    eigs = np.sort(eigs.real)[::-1]
    egv_max = max( eigs )
    
    return max(0.0,2*sqrt(egv_max)-sqrt(eigs[0])-sqrt(eigs[1])-sqrt(eigs[2])-sqrt(eigs[3]))
	