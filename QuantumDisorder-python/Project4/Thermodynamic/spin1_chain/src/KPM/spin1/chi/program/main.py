#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import math
import numpy as np
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg as splin
import scipy.fftpack as fftp

import gen_operator
import gen_hamiltonian
import time
import random
import gen_diagonal_ME as GME



#Rescaling the Hamiltonian to the interval [-1 ,1]
def rescale(H):
    dim=H.shape[0]
    lmax = float(splin.eigsh(H, k=1, which='LA',return_eigenvectors=False , ncv=25, maxiter=1E4))
    lmin = float(splin.eigsh(H, k=1, which='SA',return_eigenvectors=False , ncv=25, maxiter=1E4))
    a = (lmax - lmin) / 1.99
    b = (lmax + lmin) / 2
    H_rescaled = (1/a) * (H - b * scipy.sparse.eye(dim))
    return H_rescaled, a, b


#The Jackson Kernel improvement
def Jackson_kernel(N):
    n = np.arange(N)
    return ((N-n+1)*np.cos(np.pi*n/(N+1))+np.sin(np.pi*n/(N+1))
            /np.tan(np.pi/(N+1)))/(N+1)



def get_hamiltonian_sparse(L,K,J):
    bc = 'pbc'
    Lmax = L if bc == 'pbc' else L-1
    bonds = [(i,np.mod(i+1,L)) for i in range(Lmax)]

    Jxx = [0.5 * J] * len(bonds)
    Jyy = Jxx
    Jzz = [J] * len(bonds)


    #--------------------------------
    # %% Generate operators
    op  =  gen_operator.Operator(L)

    hamiltonian = 0.

    #--------------------------------
    #%% XXZ term
    for OP,J in zip(['sx','sy','sz'],[Jxx,Jyy,Jzz]):
        op.gen_s0sxsysz(ope=OP)
        ham = gen_hamiltonian.Hamiltonian(op.s_list)
        ham.gen_nn_int(J,bonds)
        hamiltonian +=  ham.H

    #--------------------------------
    #%% field term
    hz = [0.01] * L
    op.gen_s0sxsysz(ope='sz')
    ham = gen_hamiltonian.Hamiltonian(op.s_list)
    ham.gen_onsite_field(hz)
    hamiltonian +=  ham.field

    #--------------------------------
    #%% single-site anisotropy term
    Kz = [-K] * L
    op.gen_s0sxsysz(ope='sz')
    ham = gen_hamiltonian.Hamiltonian(op.s_list)
    ham.gen_onsite_anisotropy(Kz)
    hamiltonian +=  ham.aniotropy


    return hamiltonian


def calc_thermodynamic_cv(H_rescaled , scale_fact_a , scale_fact_b, N,num_vec,K,beta):

    """ Calculate the density of states of a Hamiltonian
        Parameters :
        ---------------------
        H_rescaled : sparse matrix
                rescaled Hamiltonian of the system

        N : integer
            Number o f Chebyshev moments
        num_vec : integer
            Number of random vectors used for sampling
        K : integer
            Number of points at which the density of states is computed
            Must be larger than N

        beta:  1.0/(kB*T)


        Returns :
        -----------------

        Cv : 1D array
                Density of states at K points
        xk : 1D array
                Points at which the density of states is evaluated
    """

    #Rescale Hamiltonian
    #H_rescaled , scale_fact_a , scale_fact_b = rescale(H)

    dim = H_rescaled.shape[0]
    #Empty vector mu for constructing the moments
    mu = np.zeros(N)

    for r in range(num_vec):
        rand_vect = np.exp(1j * 2*np.pi * np.random.rand(dim))

        #Use iterative scheme to calculate the moments
        alpha = []
        alpha.append(rand_vect)
        alpha.append(H_rescaled.dot(rand_vect))

        #mu single state is the moment using only one random vector,
        #add up all the states to get the moments mu
        mu_single_state = np.zeros(N, dtype=complex)
        mu_single_state[0] = (alpha[0].T.conj()).dot(alpha[0])
        mu_single_state[1] = (alpha[0].T.conj()).dot(alpha[1])
        for n in range(1,N//2):
            alpha.append(2*H_rescaled.dot(alpha[n])-alpha[n-1])
            #Use the symmetrical relation discussed in the section calculation
            # of moments of the density of states

            mu_single_state[2*n]   = 2*(alpha[n].T.conj()).dot(alpha[n])   - mu_single_state[0]
            mu_single_state[2*n+1] = 2*(alpha[n+1].T.conj()).dot(alpha[n]) - mu_single_state[1]
        mu = mu + mu_single_state.real
    mu = mu/num_vec/dim


    mu_ext = np.zeros(K)
    #Apply the Jackson Kernel improvement
    mu_ext[0:N] = mu * Jackson_kernel(N)
    #Use the discrete cosine trans form to get back to the original density of states
    mu_T = fftp.dct(mu_ext ,type=3)
    k = np.arange(0, K)
    #Define the Chebyshev nodes xk
    xk_rescaled = np.cos(np.pi*(k+0.5)/K)

    #Final multiplication and rescaling to get the density of  states
    gk  = np.pi * np.sqrt(1.-xk_rescaled**2)
    xk  = xk_rescaled * scale_fact_a + scale_fact_b

    rho = np.divide(mu_T ,gk)/(scale_fact_a)

    Z   = np.sum(mu_T * np.exp(-beta*xk))/K
    E   = np.sum(mu_T * np.exp(-beta*xk) * xk)/K
    E2  = np.sum(mu_T * np.exp(-beta*xk) * xk**2)/K

    return xk, rho, Z, E, E2

########################################
########################################
########################################
########################################
########################################

def calc_exp_operator(A,H_rescaled , scale_fact_a , scale_fact_b, N,num_vec,K,beta):

    """ Calculate the expectation value of operator A
        Parameters :
        ---------------------
        H_rescaled : sparse matrix
            rescaled Hamiltonian of the system
        A : sparse matrix
            Operator to be exvaluated <A>
        N : integer
            Number o f Chebyshev moments
        num_vec : integer
            Number of random vectors used for sampling
        K : integer
            Number of points at which the density of states is computed
            Must be larger than N

        beta:  1.0/(kB*T)


        Returns :
        -----------------

        <A> : 1D array
                Density of states at K points

    """


    #Rescale Hamiltonian
    #H_rescaled , scale_fact_a , scale_fact_b = rescale(H_rescaled)

    dim = H_rescaled.shape[0]
    #Empty vector mu for constructing the moments
    mu = np.zeros(N)

    for r in range(num_vec):
        rand_vect = np.exp(1j * 2*np.pi * np.random.rand(dim))

        #Use iterative scheme to calculate the moments
        alpha = []
        alpha.append(rand_vect)
        alpha.append(H_rescaled.dot(rand_vect))

        #mu single state is the moment using only one random vector,
        # add up all the states to get the moments mu
        mu_single_state = np.zeros(N, dtype=complex)
        for n in range(1,N-1):
            alpha.append(2*H_rescaled.dot(alpha[n])-alpha[n-1])
        for n in range(N):
            mu_single_state[n] = alpha[0].T.conj() @ A @ alpha[n]
        mu = mu + mu_single_state.real
    mu = mu/num_vec/dim


    mu_ext = np.zeros(K)
    #Apply the Jackson Kernel improvement
    mu_ext[0:N] = mu * Jackson_kernel(N)
    #Use the discrete cosine trans form to get back to the original density of states

    k = np.arange(0, K)
    #Define the Chebyshev nodes xk
    xk_rescaled = np.cos(np.pi*(k+0.5)/K)

    #Final multiplication and rescaling to get the density of  states
    gk  = np.pi * np.sqrt(1.0 - xk_rescaled**2)
    xk  = xk_rescaled * scale_fact_a + scale_fact_b
    expectation_values = fftp.dct(mu_ext ,type=3)#/gk

    expectation_values   = np.sum(expectation_values * np.exp(-beta*xk) )/K


    return expectation_values




########################################
########################################
########################################
########################################
L=10
k=0.0
ham = get_hamiltonian_sparse(L,k,J=1.0)
H_rescaled , scale_fact_a , scale_fact_b = rescale(ham)
#dim = ham.shape[0]
#Id  = scipy.sparse.eye(dim)


#--------------------------------
#  mesurment operators
#--------------------------------
# the squared magnetization M=(\sum_i S_i^z)/L
# \chi(T) = (<M**2> - <M>**2 )/T
hz = [1.0] * L
# %% Generate operators
op  =  gen_operator.Operator(L)
op.gen_s0sxsysz(ope='sz')
mz = gen_hamiltonian.Hamiltonian(op.s_list)
mz.gen_onsite_field(hz)
Mz =  mz.field
del mz


########################################
##### define system parameters #####
# System properties ,
N       = 100                  # number of moments
num_vec = 50                   # number of random states
K       = 2 *N                 # number of points
T       = np.logspace(-1.4,0.5,50,base=10)
#T       = np.linspace(0.01,2,100) # temperature vector
beta    = 1.0/(T+1e-15)          # inverse temperature vector

#T = np.loadtxt('temp')
#beta    = 1.0/(T+1e-15)
##### finite temperature methods #####
#
# preallocate lists to store results from iterations


f1 = open(f'thermodynamic_KPM_L{L}_K{k}.txt','w')
f1.write('Temperature,  F    <E>     <E^2>    Cv    Chi    entropy\n')
f1.write('============================================\n')

for bt in beta:
    start = time.time()
    xk, rho, Z, E, E2  = calc_thermodynamic_cv(H_rescaled , scale_fact_a , scale_fact_b,N,num_vec,K,bt)
    M                  = calc_exp_operator(Mz,H_rescaled , scale_fact_a , scale_fact_b, N,num_vec,K,bt)
    M2                 = calc_exp_operator(Mz**2,H_rescaled , scale_fact_a , scale_fact_b, N,num_vec,K,bt)
    freeEnergy =  -1.0   * np.log(Z) / bt
    entropy    = bt    *  (E/Z - freeEnergy)
    Cv         = bt**2 * (E2/Z-(E/Z)**2)
    Chi        = bt    * (M2/Z-(M/Z)**2)

    f1.write('%f   %f   %f    %f    %f   %f     %f\n'%(1./bt, freeEnergy, E/Z, E2/Z, Cv, Chi, entropy))
    f1.flush()
    print(f'time elapsed {time.time()-start}, Chi : {Chi}',flush=True)
    f1.close()

    #--------------------------------
#%%  memory check

#plt.semilogx(T,Chi,'.-b')
#plt.show()
