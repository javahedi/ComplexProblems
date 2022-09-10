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
 
 

def get_hamiltonian_sparse(L,J1,J2,hz,hx):

    J1_bond = [(0,1)  , (1,2)  , (2,3)  ,
                (4,5)  , (5,6)  , (6,7)  ,
                (8,9)  , (9,10) , (10,11),
                (12,13), (13,14), (14,15),
                (0,3)  , (4,7)  , (8,11) ,
                (12,15), (0,12) , (1,13) ,
                (2,14) , (3,15)]

    J_list = []
    for b in J1_bond:
        J_list.append(0.25 * J1 )


               
    # %% Generate operators
    op  =  gen_operator.Operator(L)

    hamiltonian = 0
    # ZZ
    for OP,J in zip(['sx','sx','sz'],[J_list,J_list,J_list]):
        op.gen_s0sxsysz(ope=OP)
        ham = gen_hamiltonian.Hamiltonian(op.s_list)
        ham.gen_nn_int(J,J1_bond)
        hamiltonian +=  ham.H
            
        # create magnetization-squared operator
        #ham.gen_magnetization_squared()
        #M2 = ham.MM
        
    # field in z direction
    op.gen_s0sxsysz(ope='sz')
    ham = gen_hamiltonian.Hamiltonian(op.s_list)
    ham.gen_onsite_field(L*[hz])
    hamiltonian +=  ham.h
    
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

    """ Calculate the density of states of a Hamiltonian
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
    #H_rescaled , scale_fact_a , scale_fact_b = rescale(H)
    H_rescaled = A @ H_rescaled
    
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
    gk  = np.pi * np.sqrt(1.0 - xk_rescaled**2)
    xk  = xk_rescaled * scale_fact_a + scale_fact_b

    expA   = np.sum(mu_T * np.exp(-beta*xk))/K

    return expA




########################################
########################################
########################################
########################################
L=16
ham = get_hamiltonian_sparse(L,J1=1.0,J2=0.0,hx=0.0,hz=0.01)

H_rescaled , scale_fact_a , scale_fact_b = rescale(ham)

dim = ham.shape[0]
Id  = scipy.sparse.eye(dim)

########################################
##### define system parameters #####
# System properties ,
N       = 100                   # number of moments
num_vec = 200                    # number of random states
K       = 2 *N                  # number of points
T       = np.logspace(-2,1,100,base=10)
#T       = np.linspace(0.01,2,100) # temperature vector
beta    = 1.0/(T+1e-15)          # inverse temperature vector


##### finite temperature methods #####
#
# preallocate lists to store results from iterations

Cv = []
f1 = open(f'thermodynamic_KPM_L{L}.txt','w')
f1.write('Temperature,  F    <E>     <E^2>    Cv    entropy\n')
f1.write('============================================\n')

cte =0
for bt in beta:
    start = time.time()
    xk, rho, Z, E, E2  = calc_thermodynamic_cv(H_rescaled , scale_fact_a , scale_fact_b,N,num_vec,K,bt)
    freeEnergy =  -1.0 * np.log(Z) / bt
    entropy = bt *  (E - freeEnergy)
    Cv.append(bt**2 * (E2/Z-(E/Z)**2))
    
    #Z  = calc_exp_operator(Id, H_rescaled , a , b , N,num_vec,K,bt)
    f1.write('%f   %f   %f    %f    %f   %f\n'%(1./bt, freeEnergy, E, E2, Cv[-1], entropy))
    f1.flush()
    cte+=1
    print(f'loop {cte} & time elapsed {time.time()-start}',flush=True)
f1.close()

    #--------------------------------
#%%  memory check
    
plt.semilogx(T,Cv,'.-b')
plt.show()

