#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 18:57:52 2020

@author: javadvahedi
"""

import random
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True
plt.rcParams['font.family'] = 'Times New Roman'


fig, (ax1, ax2) = plt.subplots(2, 1,figsize=(5,6))
ax1.tick_params(direction='in',axis='both', color='k',labelsize=14,
               left=True, top=True, right=True, bottom=True)

ax2.tick_params(direction='in',axis='both', color='k',labelsize=14,
               left=True, top=True, right=True, bottom=True)




L=200
N=20
J0=1
alpha=0.1

xx_bonds = []
for i in range(N):
    for j in range(i+1,N):  # N-1  for PBC
        xx_bonds.append((i,j))




dist=[]
J_dist=[]
for i in range(500):
    site_index = random.sample(range(L), N)
    dist = dist + site_index
    
    for bond in xx_bonds:           
            n,m=bond[0],bond[1]
            J_dist.append(J0 * abs(site_index[n]-site_index[m])**-alpha)
    
    
    





# the histogram of the data
n, bins, patches = ax1.hist(dist, 100, density=False, facecolor='g', alpha=0.75)
ax1.axes.set_xlabel(r'L',fontsize=20)
ax1.axes.set_ylabel(r'Distribution with n=0.1',fontsize=16)
ax1.axes.set_title(r'Histogram of 500 realization',fontsize=20)
ax1.axes.set_xlim(0, L)
#ax1.axes.set_ylim(0, 0.03)


n, bins, patches = ax2.hist(J_dist, 80, density=False, 
                            facecolor='b', alpha=0.75,label=r'$\alpha='+str(alpha)+'$')
ax2.axes.set_xlabel(r'J',fontsize=20)
ax2.axes.set_ylabel(r'Distribution ',fontsize=18)
ax2.axes.set_title(r'Histogram of 500 realization',fontsize=20)
ax2.axes.set_xlim(0, 1.0)
ax2.legend(loc='best',fontsize=15,frameon=False)

    
plt.tight_layout()    
plt.savefig('Distribution_alpha_01.pdf')
plt.show()
