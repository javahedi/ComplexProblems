#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 22:40:46 2018

@author: Javad Vahedi
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from numpy import pi

import matplotlib as mpl
mpl.rcParams['text.usetex']        = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.family']       = 'sans'
#mpl.rcParams['font.family']        = 'DejaVu Sans'
mpl.rcParams['lines.linewidth']    = 2
mpl.rcParams['lines.markersize']   = 6
mpl.rcParams['font.size']          = 9 
mpl.rcParams['legend.fontsize']    = 14
mpl.rcParams['axes.titlesize']     = 14
mpl.rcParams['axes.labelsize']     = 14
mpl.rcParams['xtick.major.size']   = 3
mpl.rcParams['ytick.major.size']   = 3
mpl.rcParams['xtick.major.width']  = 1
mpl.rcParams['ytick.major.width']  = 1
mpl.rcParams['xtick.direction']    = 'in'
mpl.rcParams['ytick.direction']    = 'in'
mpl.rcParams['xtick.bottom']       = 'on'
mpl.rcParams['xtick.top']          = 'on'
mpl.rcParams['ytick.left']         = 'on'
mpl.rcParams['ytick.right']        = 'on'
mpl.rcParams['figure.titlesize']   = 14


#  band
Band_ez00 = np.loadtxt('Band_ez_0.0.txt', unpack=True)
Band_ez01 = np.loadtxt('Band_ez_0.1.txt', unpack=True)
Band_ez02 = np.loadtxt('Band_ez_0.2.txt', unpack=True)

t1 = -1.22

Band_ez00/=abs(t1)
Band_ez01/=abs(t1)
Band_ez02/=abs(t1)
##############################################
##############################################
##############################################
Nx = 6      #Number of unit cell in x-direction 
Ny = 48     #Number of unit cell in y-direction 
NN = 4*Nx   # each unit cell has 4 atom
Nk=201      # number of k-point
b = 3.32
kx=np.linspace(-pi/b,pi/b,Nk)
xx = np.zeros(Nk)

my_xticks = ['$-1$', '0', '$1$']


f = plt.figure(figsize=(6,5))

ax = f.add_subplot(1,3,1)
for i in range(NN):
    if i==int(NN/2) or i==int(NN/2)-1:
        plt.plot(kx,Band_ez00[i,:],'-',color='gold')
    else:
        plt.plot(kx,Band_ez00[i,:],'-k')
plt.xticks([-pi/b,0/b,pi/b], my_xticks, fontsize = 12) 
plt.ylim([-4,4])     
plt.xlim([-pi/b,pi/b])     
plt.plot(kx,xx,'--r',linewidth=0.5)      
ax.set_xlabel(r'$kb/\pi$',fontsize=14)
#ax.axes.xaxis.set_ticklabels([])
ax.set_ylabel(r'$E-E_F(|t_1|)$',fontsize=12)
ax.tick_params(direction='in',axis='both',
               left='on', top='on', right='on', bottom='on')
plt.text(-0.5, 4.08, r'$\epsilon_z=0\%$',  color='k', fontsize=15)


ax = f.add_subplot(1,3,2)
for i in range(NN):
    if i==int(NN/2) or i==int(NN/2)-1:
        plt.plot(kx,Band_ez01[i,:],'-',color='gold')
    else:
        plt.plot(kx,Band_ez01[i,:],'-k')
plt.xticks([-pi/b,0/b,pi/b], my_xticks, fontsize = 12) 
plt.ylim([-4,4]) 
plt.xlim([-pi/b,pi/b])     
plt.plot(kx,xx,'--r',linewidth=0.5)           
ax.set_xlabel(r'$kb/\pi$',fontsize=14)
#ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both',
               left='on', top='on', right='on', bottom='on')
plt.text(-0.5, 4.08, r'$\epsilon_z=10\%$',  color='k', fontsize=15)

ax = f.add_subplot(1,3,3)
for i in range(NN):
    if i==int(NN/2) or i==int(NN/2)-1:
        plt.plot(kx,Band_ez02[i,:],'-',color='gold')
    else:
        plt.plot(kx,Band_ez02[i,:],'-k')
plt.xticks([-pi/b,0/b,pi/b], my_xticks, fontsize = 12)  
plt.ylim([-4,4])     
plt.xlim([-pi/b,pi/b])     
plt.plot(kx,xx,'--r',linewidth=0.5) 
ax.set_xlabel(r'$kb/\pi$',fontsize=14)
#ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both',
               left='on', top='on', right='on', bottom='on')

plt.text(-0.5, 4.08, r'$\epsilon_z=20\%$',  color='k', fontsize=15)





plt.tight_layout()
plt.savefig('fig6.pdf')
plt.show(f)
