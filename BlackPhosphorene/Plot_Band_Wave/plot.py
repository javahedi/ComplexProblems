#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 22:40:46 2018

@author: Javad Vahedi
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import matplotlib.cm as cm
from numpy import pi



#  down
Wave_d_U00 = np.loadtxt('Wave_down_U_0.0.txt', unpack=True)
Wave_d_U01 = np.loadtxt('Wave_down_U_0.1.txt', unpack=True)
Wave_d_U02 = np.loadtxt('Wave_down_U_0.2.txt', unpack=True)
Wave_d_U03 = np.loadtxt('Wave_down_U_0.3.txt', unpack=True)
Wave_d_U04 = np.loadtxt('Wave_down_U_0.4.txt', unpack=True)
Wave_d_U05 = np.loadtxt('Wave_down_U_0.5.txt', unpack=True)
Wave_d_U06 = np.loadtxt('Wave_down_U_0.6.txt', unpack=True)
Wave_d_U07 = np.loadtxt('Wave_down_U_0.7.txt', unpack=True)
Wave_d_U08 = np.loadtxt('Wave_down_U_0.8.txt', unpack=True)


#  up
Wave_u_U00 = np.loadtxt('Wave_up_U_0.0.txt', unpack=True)
Wave_u_U01 = np.loadtxt('Wave_up_U_0.1.txt', unpack=True)
Wave_u_U02 = np.loadtxt('Wave_up_U_0.2.txt', unpack=True)
Wave_u_U03 = np.loadtxt('Wave_up_U_0.3.txt', unpack=True)
Wave_u_U04 = np.loadtxt('Wave_up_U_0.4.txt', unpack=True)
Wave_u_U05 = np.loadtxt('Wave_up_U_0.5.txt', unpack=True)
Wave_u_U06 = np.loadtxt('Wave_up_U_0.6.txt', unpack=True)
Wave_u_U07 = np.loadtxt('Wave_up_U_0.7.txt', unpack=True)
Wave_u_U08 = np.loadtxt('Wave_up_U_0.8.txt', unpack=True)

#  band
Band_U00 = np.loadtxt('Band_U_0.0.txt', unpack=True)
Band_U01 = np.loadtxt('Band_U_0.1.txt', unpack=True)
Band_U02 = np.loadtxt('Band_U_0.2.txt', unpack=True)
Band_U03 = np.loadtxt('Band_U_0.3.txt', unpack=True)
Band_U04 = np.loadtxt('Band_U_0.4.txt', unpack=True)
Band_U04 = np.loadtxt('Band_U_0.4.txt', unpack=True)
Band_U05 = np.loadtxt('Band_U_0.5.txt', unpack=True)
Band_U0585 = np.loadtxt('Band_U_0.585.txt', unpack=True)
Band_U06 = np.loadtxt('Band_U_0.6.txt', unpack=True)
Band_U07 = np.loadtxt('Band_U_0.7.txt', unpack=True)
Band_U08 = np.loadtxt('Band_U_0.8.txt', unpack=True)

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
yy = np.zeros(Nk)

my_xticks = ['$-\\pi$', '0', '$\\pi$']


f = plt.figure(figsize=(10,10))

ax = f.add_subplot(3,5,1)
for i in range(NN):
    if min(abs(Band_U00[i,:]))> 0.4:
        plt.plot(kx,Band_U00[i,:],'-b')
    else:
        plt.plot(kx,Band_U00[i,:],'-r')
#plt.xticks([-pi/b,0/b,pi/b], my_xticks, fontsize = 12)       
#ax.set_xlabel('kb',fontsize=12)
ax.axes.xaxis.set_ticklabels([])
ax.set_ylabel(r'$E(k)-E_F[eV]$',fontsize=15)
ax.tick_params(direction='in',axis='both',
               left='on', top='on', right='on', bottom='on')

plt.xlim([-pi/b,pi/b])     
#plt.plot(kx,xx,'-k',linewidth=0.5) 

axx = plt.axes([.14, 0.86, .12, .14],)
for i in range(NN):
    if min(abs(Band_U00[i,:]))> 0.4:
        plt.plot(kx,Band_U00[i,:],'-b')
    else:
        plt.plot(kx,Band_U00[i,:],'-r')
        
axx.axes.yaxis.set_ticklabels([])
axx.axes.xaxis.set_ticklabels([])
axx.tick_params(direction='in',axis='both',
               left='on', top='on', right='on', bottom='on')
axx.set_ylim([-0.4,0.2])
plt.xlim([-pi/b,pi/b])     
plt.plot(kx,xx,'-k',linewidth=0.5) 




ax = f.add_subplot(3,5,2)
for i in range(NN):
    if min(abs(Band_U02[i,:]))> 0.4:
        plt.plot(kx,Band_U02[i,:],'-b')
    else:
        plt.plot(kx,Band_U02[i,:],'-r')
#plt.xticks([-pi/b,0/b,pi/b], my_xticks, fontsize = 12)       
#ax.set_xlabel('kb',fontsize=12)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both',
               left='on', top='on', right='on', bottom='on')
plt.xlim([-pi/b,pi/b])     
#plt.plot(kx,xx,'-k',linewidth=0.5) 




axx = plt.axes([.32, .86, .12, .14],)
for i in range(NN):
    if min(abs(Band_U02[i,:]))> 0.4:
        plt.plot(kx,Band_U02[i,:],'-b')
    else:
        plt.plot(kx,Band_U02[i,:],'-r')
        
axx.axes.xaxis.set_ticklabels([])        
axx.axes.yaxis.set_ticklabels([])
axx.tick_params(direction='in',axis='both',
               left='on', top='on', right='on', bottom='on')
axx.set_ylim([-0.4,0.2])
plt.xlim([-pi/b,pi/b])     
plt.plot(kx,xx,'-k',linewidth=0.5) 




ax = f.add_subplot(3,5,3)
for i in range(NN):
    if min(abs(Band_U04[i,:]))> 0.4:
        plt.plot(kx,Band_U04[i,:],'-b')
    else:
        plt.plot(kx,Band_U04[i,:],'-r')

#plt.xticks([-pi/b,0/b,pi/b], my_xticks, fontsize = 12)       
#ax.set_xlabel('kb',fontsize=12)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both',
               left='on', top='on', right='on', bottom='on')
plt.xlim([-pi/b,pi/b])     
#plt.plot(kx,xx,'-k',linewidth=0.5) 


axx = plt.axes([.49, .86, .12, .14],)
for i in range(NN):
    if min(abs(Band_U04[i,:]))> 0.4:
        plt.plot(kx,Band_U04[i,:],'-b')
    else:
        plt.plot(kx,Band_U04[i,:],'-r')
        
axx.axes.xaxis.set_ticklabels([])        
axx.axes.yaxis.set_ticklabels([])
axx.tick_params(direction='in',axis='both',
               left='on', top='on', right='on', bottom='on')
axx.set_ylim([-0.16,0.5])
plt.xlim([-pi/b,pi/b])     
plt.plot(kx,xx,'-k',linewidth=0.5)



ax = f.add_subplot(3,5,4)
for i in range(NN):
    if min(abs(Band_U06[i,:]))> 0.4:
        plt.plot(kx,Band_U06[i,:],'-b')
    else:
        plt.plot(kx,Band_U06[i,:],'-r')

#plt.xticks([-pi/b,0/b,pi/b], my_xticks, fontsize = 12)       
#ax.set_xlabel('kb',fontsize=12)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both',
               left='on', top='on', right='on', bottom='on')
plt.xlim([-pi/b,pi/b])     
#plt.plot(kx,xx,'-k',linewidth=0.5) 

axx = plt.axes([.66, .86, .12, .14],)
for i in range(NN):
    if min(abs(Band_U0585[i,:]))> 0.4:
        plt.plot(kx,Band_U0585[i,:],'-b')
    else:
        plt.plot(kx,Band_U0585[i,:],'-r')
        
axx.axes.xaxis.set_ticklabels([])        
axx.axes.yaxis.set_ticklabels([])
axx.tick_params(direction='in',axis='both',
               left='on', top='on', right='on', bottom='on')
axx.set_ylim([-0.16,0.5])
plt.xlim([-pi/b,pi/b])     
plt.plot(kx,xx,'-k',linewidth=0.5)

rectangle = patches.Rectangle((-pi/b,0),         # (x,y)
                       2*pi/b,          # width
                       min(abs(Band_U0585[int(NN/2),:])),          # height
                       #fill=False,
                       edgecolor='None',
                       facecolor='yellow',
                       alpha = 1.0)
axx.add_patch(rectangle)



ax = f.add_subplot(3,5,5)
for i in range(NN):
    if min(abs(Band_U08[i,:]))> 0.4:
        plt.plot(kx,Band_U08[i,:],'-b')
    else:
        plt.plot(kx,Band_U08[i,:],'-r')

#plt.xticks([-pi/b,0/b,pi/b], my_xticks, fontsize = 12)       
#ax.set_xlabel('kb',fontsize=12)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both',
               left='on', top='on', right='on', bottom='on')
plt.xlim([-pi/b,pi/b])     
#plt.plot(kx,xx,'-k',linewidth=0.5) 


axx = plt.axes([.845, .86, .12, .14],)
for i in range(NN):
    if min(abs(Band_U08[i,:]))> 0.4:
        plt.plot(kx,Band_U08[i,:],'-b')
    else:
        plt.plot(kx,Band_U08[i,:],'-r')
      
        
        
axx.axes.xaxis.set_ticklabels([])        
axx.axes.yaxis.set_ticklabels([])
axx.tick_params(direction='in',axis='both',
               left='on', top='on', right='on', bottom='on')
axx.set_ylim([-0.16,0.5])
plt.xlim([-pi/b,pi/b])     
plt.plot(kx,xx,'-k',linewidth=0.5)

#plt.tight_layout()
#plt.savefig('Fig5a.pdf')
#plt.show(f1)
###########
###########
###########
###########
###########
#f2 = plt.figure(figsize=(8,6))

ax = f.add_subplot(3,5,6)
plt.imshow(Wave_u_U00.T, interpolation='bilinear',
           cmap=cm.gnuplot,origin='lower',
           extent=[-5, 5, 0, NN])

my_yticks = ['0', '12', '24']
plt.yticks([0,int(NN/2),NN], my_yticks, fontsize = 10)
plt.ylabel(r'$x(\AA)$',fontsize=18)
#ax.set_ylabel(r'$y(\AA)$',fontsize=15)
ax.axes.xaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both', color='w',
               left='on', top='on', right='on', bottom='on')
plt.text(-2, 11, r'$\|\uparrow\rangle$',  color='w', fontsize=20)
ax.set_aspect('auto')



ax = f.add_subplot(3,5,7)
plt.imshow(Wave_u_U02.T, interpolation='bilinear',
           cmap=cm.gnuplot,origin='lower',
           extent=[-5, 5, 0, NN])

plt.yticks([0,int(NN/2),NN], ['','',''], fontsize = 10)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both', color='w',
               left='on', top='on', right='on', bottom='on')
plt.text(-2, 11, r'$\|\uparrow\rangle$',  color='w', fontsize=20)
ax.set_aspect('auto')




ax = f.add_subplot(3,5,8)
plt.imshow(Wave_u_U04.T, interpolation='bilinear',
           cmap=cm.gnuplot,origin='lower',
           extent=[-5, 5, 0, NN])

plt.yticks([0,int(NN/2),NN], ['','',''], fontsize = 10)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both', color='w',
               left='on', top='on', right='on', bottom='on')
plt.text(-2, 11, r'$\|\uparrow\rangle$',  color='w', fontsize=20)
ax.set_aspect('auto')




ax = f.add_subplot(3,5,9)
plt.imshow(Wave_u_U06.T, interpolation='bilinear',
           cmap=cm.gnuplot,origin='lower',
           extent=[-5, 5, 0, NN])

plt.yticks([0,int(NN/2),NN], ['','',''], fontsize = 10)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both', color='w',
               left='on', top='on', right='on', bottom='on')
plt.text(-2, 11, r'$\|\uparrow\rangle$',  color='w', fontsize=20)
ax.set_aspect('auto')




ax = f.add_subplot(3,5,10)
cax =plt.imshow(Wave_u_U08.T, interpolation='bilinear',
           cmap=cm.gnuplot,origin='lower',
           extent=[-5, 5, 0, NN])
cbar = f.colorbar(cax, ticks=[0, 0.9], orientation='vertical',)
cbar.ax.set_xticklabels(['0', '1'])  # horizontal colorbar
ax.set_aspect('auto')
plt.yticks([0,int(NN/2),NN], ['','',''], fontsize = 10)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both', color='w',
               left='on', top='on', right='on', bottom='on')
plt.text(-2, 11, r'$\|\uparrow\rangle$',  color='w', fontsize=20)

############
############
############
############
ax = f.add_subplot(3,5,11)
plt.imshow(Wave_d_U00.T, interpolation='bilinear',
           cmap=cm.gnuplot,origin='lower',
           extent=[-5, 5, 0, NN])

my_yticks = ['0', '12', '24']
plt.yticks([0,int(NN/2),NN], my_yticks, fontsize = 10)
plt.ylabel(r'$x(\AA)$',fontsize=18)
plt.xticks([-5,0,5], my_xticks, fontsize = 15) 
ax.set_xlabel('kb',fontsize=15)
ax.tick_params(direction='in',axis='both', color='w',
               left='on', top='on', right='on', bottom='on')
plt.text(-2, 11, r'$\|\downarrow\rangle$',  color='w', fontsize=20)
ax.set_aspect('auto')

ax = f.add_subplot(3,5,12)
plt.imshow(Wave_d_U02.T, interpolation='bilinear',
           cmap=cm.gnuplot,origin='lower',
           extent=[-5, 5, 0, NN])

plt.yticks([0,int(NN/2),NN], ['','',''], fontsize = 10)
plt.xticks([-5,0,5], my_xticks, fontsize = 15) 
ax.set_xlabel('kb',fontsize=15)
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both', color='w',
               left='on', top='on', right='on', bottom='on')
plt.text(-2, 11, r'$\|\downarrow\rangle$',  color='w', fontsize=20)
ax.set_aspect('auto')


ax = f.add_subplot(3,5,13)
plt.imshow(Wave_d_U04.T, interpolation='bilinear',
           cmap=cm.gnuplot,origin='lower',
           extent=[-5, 5, 0, NN])

plt.yticks([0,int(NN/2),NN], ['','',''], fontsize = 10)
plt.xticks([-5,0,5], my_xticks, fontsize = 15) 
ax.set_xlabel('kb',fontsize=15)
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both', color='w',
               left='on', top='on', right='on', bottom='on')
plt.text(-2, 11, r'$\|\downarrow\rangle$',  color='w', fontsize=20)
ax.set_aspect('auto')


ax = f.add_subplot(3,5,14)
plt.imshow(Wave_d_U06.T, interpolation='bilinear',
           cmap=cm.gnuplot,origin='lower',
           extent=[-5, 5, 0, NN])

plt.yticks([0,int(NN/2),NN], ['','',''], fontsize = 10)
plt.xticks([-5,0,5], my_xticks, fontsize = 15) 
ax.set_xlabel('kb',fontsize=15)
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both', color='w',
               left='on', top='on', right='on', bottom='on')
plt.text(-2, 11, r'$\|\downarrow\rangle$',  color='w', fontsize=20)
ax.set_aspect('auto')



ax = f.add_subplot(3,5,15)
cax = ax.imshow(Wave_d_U08.T, interpolation='bilinear',
           cmap=cm.gnuplot,origin='lower',
           extent=[-5, 5, 0, NN])
cbar = f.colorbar(cax, ticks=[0, 0.9], orientation='vertical')
cbar.ax.set_xticklabels(['0', '1'])  # horizontal colorbar
ax.set_aspect('auto')
plt.yticks([0,int(NN/2),NN], ['','',''], fontsize = 10)
plt.xticks([-5,0,5], my_xticks, fontsize = 15) 
ax.set_xlabel('kb',fontsize=15)
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(direction='in',axis='both', color='w',
               left='on', top='on', right='on', bottom='on')
plt.text(-2, 11, r'$\|\downarrow\rangle$',  color='w', fontsize=20)




plt.tight_layout()
plt.savefig('Fig5.pdf')
plt.show(f)
