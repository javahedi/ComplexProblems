#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 12:23:49 2020

@author: javadvahedi
"""

import numpy as np
import scipy
import scipy.sparse as sparse

def gen_diagprojector(symvec, symval):
    ind = np.where(symvec==float(symval))
    dim = np.size(ind)
    P = sparse.lil_matrix((dim,len(symvec)))
    for j in range(dim):
        P[j,ind[0][j]] = 1.0
    return P