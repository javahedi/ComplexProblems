#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 20:52:44 2019

@author: https://github.com/cvxgrp/cvxpy/issues/563
"""

import numpy as np


def get_partial_trace(rho, dims, axis=0):
    """
    Takes partial trace over the subsystem defined by 'axis'
    rho: a matrix
    dims: a list containing the dimension of each subsystem
    axis: the index of the subsytem to be traced out
    (We assume that each subsystem is square)
    """
    dims_ = np.array(dims)
    # Reshape the matrix into a tensor with the following shape:
    # [dim_0, dim_1, ..., dim_n, dim_0, dim_1, ..., dim_n]
    # Each subsystem gets one index for its row and another one for its column
    reshaped_rho = rho.reshape(np.concatenate((dims_, dims_), axis=None))

    # Move the subsystems to be traced towards the end
    reshaped_rho = np.moveaxis(reshaped_rho, axis, -1)
    reshaped_rho = np.moveaxis(reshaped_rho, len(dims)+axis-1, -1)

    # Trace over the very last row and column indices
    traced_out_rho = np.trace(reshaped_rho, axis1=-2, axis2=-1)

    # traced_out_rho is still in the shape of a tensor
    # Reshape back to a matrix
    dims_untraced = np.delete(dims_, axis)
    rho_dim = np.prod(dims_untraced)
    return traced_out_rho.reshape([rho_dim, rho_dim])


"""
Test out the partial_trace numpy module by creating a matrix
rho_ABC = rho_A \otimes rho_B \otimes rho_C
Each rho_i is normalized, i.e. Tr(rho_i) = 1
"""

'''
# Generate the results we want
rho_A = np.random.rand(4, 4) + 1j*np.random.rand(4, 4)
rho_A /= np.trace(rho_A)
rho_B = np.random.rand(3, 3) + 1j*np.random.rand(3, 3)
rho_B /= np.trace(rho_B)
rho_C = np.random.rand(2, 2) + 1j*np.random.rand(2, 2)
rho_C /= np.trace(rho_C)

rho_AB = np.kron(rho_A, rho_B)
rho_BC = np.kron(rho_B, rho_C)
rho_AC = np.kron(rho_A, rho_C)
rho_ABC = np.kron(rho_AB, rho_C)

# Try to get the results by doing partial_trace's from rho_ABC
rho_AB_test = partial_trace(rho_ABC, [4, 3, 2], axis=2)
rho_AC_test = partial_trace(rho_ABC, [4, 3, 2], axis=1)
rho_A_test = partial_trace(rho_AB_test, [4, 3], axis=1)
rho_B_test = partial_trace(rho_AB_test, [4, 3], axis=0)
rho_C_test = partial_trace(rho_AC_test, [4, 2], axis=0)

# See if the outputs of partial_trace are correct
print("rho_AB test correct? ", np.allclose(rho_AB_test, rho_AB))
print("rho_AC test correct? ", np.allclose(rho_AC_test, rho_AC))
print("rho_A test correct? ", np.allclose(rho_A_test, rho_A))
print("rho_B test correct? ", np.allclose(rho_B_test, rho_B))
print("rho_C test correct? ", np.allclose(rho_C_test, rho_C))
'''