#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 15:55:52 2019

@author: Javad Vahedi & Stefan Kettemann
"""

import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
from scipy.sparse import linalg
import matplotlib.pyplot as plt
import Hamiltonian
import Operator
import partialtrace
import math
#import itertools


#plt.rcParams['text.usetex'] = True
#plt.rcParams['text.latex.unicode'] = True
#plt.rcParams["font.family"] = "Times New Roman"



# Set parameter
N  = 10#np.loadtxt('N.txt')    # number of spin
L  = 100#np.loadtxt('L.txt')   # latacice size
n  = N/L  # density
J0 = 1
hz = J0*0
alpha = 1# np.loadtxt('alpha.txt')
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

print("ground state energy (per site):", energy/N, flush=True)

del(hamiltonian)
del(energy)

for i in range(len(space_index)):
    f1.write("%i\n"  %(space_index[i]))
    f1.flush()
f1.close()


groundstate = groundstate[:,0]
sxy_corrs = []
for n in range(N):
        for m in range(n+1,N):  # N-1  for PBC
            l = abs(space_index[n] -space_index[m])
            
            sxy_rows, sxy_cols, sxy_data = Operator.get_spinspin_sparse(L,N,n,m,hz,False)
            sxy_matrix = sp.sparse.csr_matrix((sxy_data, (sxy_rows, sxy_cols)))
            sxy_corrs.append(np.dot(groundstate, sxy_matrix.dot(groundstate)))
            
            f2.write("%i  %i  %i  %i  %f\n"  %(n+1,m+1, abs(n-m), l ,sxy_corrs[-1]))
            f2.flush()
f2.close()


del(sxy_matrix)


n,m,d,l,sxx=np.loadtxt('Sxx_L'+str(L)+'_N_'+str(N)+'_a_'+str(alpha)+'.txt',unpack=True)

#for counter, value in enumerate(d.astype(int)): 
    #print((counter, value),' ', int(l[counter]),'  ' ,abs(sxx[counter]))

#print('######################')
#print('######################')


idx   = np.argsort(l)
l = np.array(l.astype(int))[idx]
sxx = np.array(abs(sxx))[idx]


from collections import Counter
counter = Counter(l)
result = [(i,j) for i, j in counter.items() if j > 1]

sxx_new = []
l_new = []
ll_index = []
for i,j in result:
    ll = np.where(l==i)[0]
    s = 0
    for n in range(j):
        s += sxx[ll[n]]
        ll_index.append(ll[n])
    sxx_new.append(s/j)
    l_new.append(l[ll[n]])


for i in range(len(l)):
    if i not in ll_index:
        sxx_new.append(sxx[i])
        l_new.append(l[i])
            
   

l_new = np.asanyarray(l_new)
sxx_new = np.asanyarray(sxx_new)


idx   = np.argsort(l_new)
l_new = np.array(l_new)[idx]
sxx_new = np.array(sxx_new[idx])


'''
v = groundstate.reshape([2**N, -1])
rho = v@v.T
del(v)
i=0
Svn = []
while i<N:
    dim  =  (N-i)*[2]
    rho  =  partialtrace.get_partial_trace(rho,dim,axis=0)
    eigs =  sp.linalg.eigh(rho, eigvals_only=True)
    svn  = -np.sum(np.real(eigs*np.log(eigs)))
    if math.isnan(svn):
        Svn.append(0)
    else:
        Svn.append(-np.sum(np.real(eigs*np.log2(eigs))))
    i +=1
    
    f3.write("%i  %f\n"  %(i, Svn[-1]))
    f3.flush()
f3.close()

plt.figure()
plt.plot(Svn)
'''

fig = plt.figure(figsize=(4,4))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],facecolor='w')
ax.tick_params(direction='in',axis='both', color='k',labelsize=14,
               left='on', top='on', right='on', bottom='on')

ax.loglog(sxx, "o-b")
ax.loglog(sxx_new, "^-r")
x=np.arange(1,len(sxx))
#ax.loglog(x,np.exp(-x),'-')
ax.loglog(x,1/(x**alpha),'-c')

ax.axes.set_xlabel(r"$|n-m|$",fontsize=18)
ax.axes.set_ylabel(r"$\langle S^y_nS^y_m\rangle$",fontsize=18)
plt.show()




