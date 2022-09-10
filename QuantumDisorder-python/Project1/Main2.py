#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 15:55:52 2019

@author: Javad Vahedi & Stefan Kettemann


python -m numpy.f2py -c partial_trace.f90 partial_trace

or 

f2py3 -c -m  partial_trace partial_trace.f90 


import partial_trace
print(partial_trace.__doc__)
"""



import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
from scipy.sparse import linalg
import matplotlib.pyplot as plt
import Hamiltonian
import Operator
#import partialtrace
import math
import ptrace
import partial_trace
import Concurence
import entanglement
import time

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True
plt.rcParams['font.family'] = 'Times New Roman'




start = time.time()

# Set parameter
N  = 12      #np.loadtxt('N.txt')    # number of spin
L  = 120     #np.loadtxt('L.txt')    # lattice size
n  = N/L     # density
J0 = 1
hz = J0*0
alpha = 1    # np.loadtxt('alpha.txt')
n_lowest_eigenvalues = 1

f1 = open('Space_Index_L'+str(L)+'_N_'+str(N)+'_a_'+str(alpha)+'.txt', 'w')
f2 = open('Sxx_L'+str(L)+'_N_'+str(N)+'_a_'+str(alpha)+'.txt', 'w')
f3 = open('Svn_L'+str(L)+'_N_'+str(N)+'_a_'+str(alpha)+'.txt', 'w')



# Step 1: evaluate GS 
rows, cols, data, space_index = Hamiltonian.get_hamiltonian_sparse(L, N, J0, hz, alpha)
hamiltonian = sp.sparse.csr_matrix((data, (rows, cols)))
energy, groundstate = sp.sparse.linalg.eigsh(hamiltonian, 
                                             k=n_lowest_eigenvalues,
                                             which='SA',
                                             return_eigenvectors=True,            
                                             maxiter=1000)

print("ground state energy (per site):", energy[0]/N, flush=True)
groundstate = groundstate[:,0]



plt.figure()
plt.spy(hamiltonian,markersize=0.5,marker='.')
plt.show()


bin = np.linspace(-60,0, num=80)
plt.figure()
n, bins, patches = plt.hist(np.log(abs(groundstate)), bin, normed=True, 
                            histtype='step',lw=2)
plt.xlabel(r"$\log{(||\phi_0\rangle|)}$",fontsize=18)
plt.show()



plt.figure()
plt.plot(abs(groundstate)**2)
plt.show()

groundstate[abs(groundstate) < 1e-6] = 0


plt.figure()
plt.plot(abs(groundstate)**2)
plt.show()


bin = np.linspace(-40,0, num=80)
plt.figure()
n, bins, patches = plt.hist(np.log(abs(groundstate)), bin, normed=True, 
                            histtype='step',lw=2)
plt.xlabel(r"$\log{(||\phi_0\rangle|)}$",fontsize=18)
plt.show()


del(hamiltonian)
del(energy)




sxy_corrs = []
for n in range(1,N):
    sxy_rows, sxy_cols, sxy_data = Operator.get_spinspin_sparse(L,N,0,n,hz,False)
    sxy_matrix = sp.sparse.csr_matrix((sxy_data, (sxy_rows, sxy_cols)))
    sxy_corrs.append(np.dot(groundstate, sxy_matrix.dot(groundstate)))
    f2.write("%i  %i  %i  %f\n"  %(1,n+1, space_index[n] ,sxy_corrs[-1]))
    f2.flush()
    
    print(sxy_corrs[-1])
f2.close()


del(sxy_matrix)


n,m,d,sxx=np.loadtxt('Sxx_L'+str(L)+'_N_'+str(N)+'_a_'+str(alpha)+'.txt',unpack=True)




####################################
Psi = groundstate.reshape([2**N, -1])
rho = Psi@Psi.T
#rho = sp.sparse.csr_matrix(Psi@Psi.T)

#plt.figure()
#plt.spy(rho,markersize=0.5,marker='.')
#plt.show()    
    

####################################
#  concurence
####################################
'''
di  = [2]*N     # Vector specifying the dimensions of the sub-systems
dr  = 2**2
Con = []
for i in range(1,N):
    ssys = [0]*N
    ssys[0]=1;ssys[i]=1     # Vector (with components equal to 0 or 1) specifying 
                            # the subsystems to be traced out.
                            # If ssys(j) = 0 the j-th subsystem is traced out. 
                            # If ssys(j) = 1 it is NOT traced out.
                            
    #   partial_trace_he(rho,di,ssys,rhor,d=shape(rho,0),nss=len(di),dr=shape(rhor,0))
    rhor = partial_trace.partial_trace_he(rho,di,ssys,dr)        
    Con.append(entanglement.concurrence_2qb(rhor))
    #Con.append(Concurence.concurrence(rhor))
    print('concurence between 1 and %i, Cons =%f'%(i+1,Con[-1]))
  
    
'''
####################################
#  Von Numann entropy
####################################
nss = N         # Number of subsystems
d   = 2**N      # Total dimension
di  = [2]*N     # Vector specifying the dimensions of the sub-systems
Svn = []
for i in range(1,N):
    
    #ssys=[0]*i+[1]*(N-i)        # Vector (with components equal to 0 or 1) specifying 
                                # the subsystems to be traced out.
                                # If ssys(j) = 0 the j-th subsystem is traced out. 
                                # If ssys(j) = 1 it is NOT traced out.
                                
    #dr  = 2**(N-i)              # Dimension of the reduced matrix
                                # (is the product of the dimensions of the 
                                # sub-systems we shall not trace out)
                                
    #rhor = partial_trace.partial_trace_he(rho,di,ssys,dr)
    
    db = 2
    da = 2**(N-i)
    rho = partial_trace.partial_trace_b_he(rho, da, db)
    #rho = ptrace.get_ptrace(rho, da, db)
       
    
    
    eigs = np.linalg.eigvalsh(rho)
    #eigs = sp.linalg.eigh(rho.toarray(),eigvals_only=True)                                  
    
    eigs[eigs <= 0] = 1e-20
    svn  = -np.sum(np.real(eigs*np.log(eigs)))
    if math.isnan(svn):
        Svn.append(0)
    else:
        Svn.append(-np.sum(np.real(eigs*np.log(eigs))))
    print('Across bond b=%i, SvN =%f'%(i,Svn[-1]))
    
    
    
    

    
   

"""
for i in range(1,N+1):
    a=[0]*i+[1]*(N-i)
    b=a[::-1]
    print(a+[1]+b)
"""
    
    
    


fig = plt.figure(figsize=(4,4))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],facecolor='w')
ax.tick_params(direction='in',axis='both', color='k',labelsize=14,
               left='on', top='on', right='on', bottom='on')

ax.loglog(abs(sxx), "o-b")
x=np.arange(1,len(sxx))
ax.loglog(x,np.exp(-x),'-r')
ax.loglog(x,1/(x**alpha),'-c')
ax.axes.set_xlim([1,15])
ax.axes.set_xlabel(r"$|n-m|$",fontsize=18)
ax.axes.set_ylabel(r"$\langle S^y_nS^y_m\rangle$",fontsize=18)
plt.show()


#print(time.time()-start)


