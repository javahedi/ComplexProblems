#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 19:41:11 2020

@author: javadvahedi
"""
import numpy as np
from numpy import exp
import matplotlib.pyplot as plt
import Structure
import Hamiltonian

import sys
sys.path.append("/Users/javadvahedi/Dropbox/TwistedBilayer/kpmpy-master/src") # path to KPM library
import kpm  # library with the kernel polynomial method

import matplotlib as mpl
mpl.rcParams['text.usetex']        = True
mpl.rcParams['text.latex.unicode'] = True
#mpl.rcParams['font.family']       = 'sans-serif'
mpl.rcParams['font.family']        = 'sans'
mpl.rcParams['lines.linewidth']    = 2
mpl.rcParams['lines.markersize']   = 6
mpl.rcParams['font.size']          = 14 
mpl.rcParams['legend.fontsize']    = 20
mpl.rcParams['axes.titlesize']     = 18
mpl.rcParams['axes.labelsize']     = 18
mpl.rcParams['xtick.major.size']   = 6
mpl.rcParams['ytick.major.size']   = 6
mpl.rcParams['xtick.major.width']  = 1
mpl.rcParams['ytick.major.width']  = 1
mpl.rcParams['xtick.direction']    = 'in'
mpl.rcParams['ytick.direction']    = 'in'
mpl.rcParams['xtick.bottom']       = True
mpl.rcParams['xtick.top']          = True
mpl.rcParams['ytick.left']         = True
mpl.rcParams['ytick.right']        = True
mpl.rcParams['figure.titlesize']   = 20
mpl.rcParams['figure.figsize']     = [5,5]





###############################################
##############################################
m,n=8,9
t0=-0.75 

filename,r1,r2,b2,b1 = Structure.UnitCell(m,n)

H00,HL1,HL2,HL1L2,HL1_L2 =\
    Hamiltonian.Hamil(t0,filename,m,n,r1,r2)
    
    


G = np.array([0.0,  0.0])
K = ((1/3)*b1+(2/3)*b2)
M =  (b1+b2)*0.5

mesh = 10
GKx = np.linspace(G[0],K[0],mesh,endpoint=True)
GKy = np.linspace(G[1],K[1],mesh,endpoint=True)
KMx = np.linspace(K[0],M[0],mesh,endpoint=True)
KMy = np.linspace(K[1],M[1],mesh,endpoint=True)


kxList = np.concatenate((GKx, KMx), axis=0)
kyList = np.concatenate((GKy, KMy), axis=0)


#plt.scatter(kxList,kyList,s=2,c='b')
#plt.quiver([0,0], [0,0], [b1[0],b2[0]], 
#           [b1[1],b2[1]], angles='xy', 
#           scale_units='xy', scale=1)
#plt.xlim(-0.2, 0.2)
#plt.ylim(-0.2, 0.2)

###############################################
###############################################
    


    
site = 223 # the entry where you want the dos, n/2 is an atom in the middle
#scale = 4.0 # this number has to be such that max(|eigenvalues|)/scale < 1
npol = 200  # number of polynomials, energy resolution goes as 1/npol
ne = 50    # number of energies to calculate (between -scale and scale)



Dim = H00.shape[0]
EnergyBand=np.zeros((len(kxList),Dim))
LDOS1 = 0
LDOS2 = 0
i=0
for kx, ky in zip(kxList,kyList):
    print(i)
    
    K = np.array([kx,ky])
    KL1,  KL2    = np.dot(K,r1),    np.dot(K,r2)
    KL1L2,KL1_L2 = np.dot(K,r1+r2), np.dot(K,r1-r2)

    ################################################
    #######            H(k)              ###########
    ################################################
     
    HH = H00+HL1.T   *exp(-1j*KL1)   +HL1   *exp(1j*KL1)+\
      	     HL2.T   *exp(-1j*KL2)   +HL2   *exp(1j*KL2)+\
             HL1L2.T *exp(-1j*KL1L2) +HL1L2 *exp(1j*KL1L2)+\
             HL1_L2.T*exp(-1j*KL1_L2)+HL1_L2*exp(1j*KL1_L2)
              
    en = np.linalg.eigvals(HH)
    idx = en.argsort()[::]   
    EnergyBand[i,:] = en[idx].real
    
    scale = max(abs(EnergyBand[i,:]))
    (x1,y1) = kpm.ldos0d(HH,i=site,scale=scale,npol=npol,ne=ne) 
    LDOS1 +=y1
    
    
    # use matrix inversion
    #x2 = np.linspace(min(EnergyBand[i,:]),max(EnergyBand[i,:]),400) # energies
    #delta = .04 # smearing
    #y2 = [((np.asmatrix(HH) - (x+1j*delta)*np.eye(Dim)).I)[site,site].imag for x in x2]
    #y2 = np.array(y2)/np.pi # pi factor
    #LDOS2 +=y2
    
    i +=1
    

LDOS2 /=i
    
    
    

plt.figure()
for i in range(Dim):    
        plt.plot(EnergyBand[:,i],'-b')

plt.ylabel(r'$E-E_F[eV]$', fontsize = 20)
my_xticks = [r'$\Gamma$', r'$K$', r'$M$']
plt.xticks([0, mesh, 2*mesh-1], my_xticks, fontsize = 15)
plt.ylim([0.1,0.4])
plt.xlim([0, 2*mesh-1])




