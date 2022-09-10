#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 15:16:02 2019

@author: javadvahedi
"""

import numpy as np
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg as splin
import scipy.fftpack as fftp
import time
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True
plt.rcParams["font.family"] = "Times New Roman"


def get_hamiltonian_sparse(L,t,e):    
    ii = np.arange(0,L-2) # indexes from 0 to n-2
    cols = np.concatenate([ii,ii+1]) # index for rows
    rows = np.concatenate([ii+1,ii]) # index for cols
    data = np.zeros(cols.shape[0]) + t # tight binding parameter    
    return sp.sparse.csc_matrix((data,(rows,cols)),shape=(L,L))
    


#Rescaling the Hamiltonian to the interval [-1 ,1]
def rescale(H):
    dim=H.shape[0]
    lmax = float(splin.eigsh(H, k=1, which='LA',return_eigenvectors=False , ncv=25))
    lmin = float(splin.eigsh(H, k=1, which='SA',return_eigenvectors=False , ncv=25))
    a = (lmax - lmin) / 1.99
    b = (lmax + lmin) / 2
    H_rescaled = (1/a) * (H - b * scipy.sparse.eye(dim))
    return H_rescaled, a, b




#The Jackson Kernel improvement
def Jackson_kernel(N):
    n = np.arange(N)
    return ((N-n+1)*np.cos(np.pi*n/(N+1))+np.sin(np.pi*n/(N+1))
            /np.tan(np.pi/(N+1)))/(N+1)
    
    
t0 = time.clock()
def calc_DoS(H,N,num_vec,K,i,ldos=True):
    
    
    """ Calculate the density of states of a Hamiltonian
        Parameters :
        ---------------------
        H : sparse matrix     
            Hamil tonian of the system
        N : integer          
            Number o f Chebyshev moments
        num_vec : integer
            Number of random vectors used for sampling
        K : integer
            Number of points at which the density of states is computed
            Must be larger than N
        i : site for Local DOS
        
        Returns :
        -----------------
    
        rho : 1D array
                Density of states at K points
        xk : 1D array
                Points at which the density of states is evaluated
    """
        
    #Rescale Hamiltonian
    H_rescaled , scale_fact_a , scale_fact_b = rescale(H)
    dim =H.shape[0]
    #Empty vector mu for constructing the moments
    mu = np.zeros(N)
    
    for r in range(num_vec):
        if ldos==True:
            rand_vect = np.exp(1j * 2*np.pi * np.random.rand(dim)) * 0.0
            rand_vect[i] = 1.0 # vector only in site i
        else: # tdos
            #Make a random vector
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
    gk = np.pi*np.sqrt(1.-xk_rescaled**2)
    xk = xk_rescaled * scale_fact_a + scale_fact_b
    rho = np.divide(mu_T,gk)/(scale_fact_a)
    return xk,rho




'''
def calc_exp_values_operator(H, A, N, num vec, K, T=0, dt=1, temp=0, integrate = True, E fermi=0):
    """
    Calculate time−dependent expectation values
    Parameters :
    −−−−−−−−−−−
    H : function or sparse matrix
        When H is a function H(t) returns the Hamiltonian at time t
    A : function or sparse matrix
        Operator of which the time dependent expectation values are computed
        WhenA is a function A(t,psi 1,psi 2) must return < psi_1 | A(t) | psi_2>,
    N : integer
        Number of Chebyshev moments.
    num vec : integer
        Number of random vectors used for sampling.
    K : integer
        Number of points for which the expectation value of A is calculated
        Must be larger than N
    T : float
        Time interval , optional , default = 0
    dt : float
        Time step , optional , default = 1
    integrate : boolean
        If integrate = True returns the integrated version of
        the expanded function using Fermi−Dirac weight temp : float
        Initial temperature of the system , optional , default = 0
    Returns :
    −−−−−−−−
    Expectation values : dense matrix
    Time−dependent expectation values , if integrate = false
    returns expectation values as a function of energy
    """
'''


       
        







# System properties ,
N       = 100   # number of moments
num_vec = 20    # number of random states
K       = 2 *N  # number of points
L       = 1000  # length of chain
t       = -1    # hopping term
e       = 0     # on-site



H = get_hamiltonian_sparse(L,t,e)
xk, rho = calc_DoS(H,N,num_vec,K,i=1,ldos=True)

t1 = time.clock()
print("Time in computation",t1-t0)

plt.figure()
plt.tick_params(direction='in',axis='both', color='k',
               left='on', top='on', right='on', bottom='on')

plt.plot(xk,rho,'.-b',lw=1,alpha=0.5)

#plt.legend(loc='best',fontsize=14)
#plt.xlim([-0.6,0.6])
#plt.ylim([1e-7,1e1])
#plt.title('Transmission N96 (meta)',fontsize=16)
#plt.savefig('Transmission2.png')
plt.show()
