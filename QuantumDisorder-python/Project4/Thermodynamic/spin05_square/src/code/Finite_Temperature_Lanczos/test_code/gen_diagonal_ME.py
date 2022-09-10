#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 20:01:08 2020

@author: javadvahedi
"""
import numpy as np


###  Assumes real eigenvectors  #### 
# (Computation is faster if you know your eigenvectors are real)
def gen_diagonal_ME(op, evecs):
    return np.sum(evecs@op@evecs,0)