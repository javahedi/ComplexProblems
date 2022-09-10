#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 21:50:58 2020

@author: javadvahedi
"""


import numpy as np
#from mayavi.mlab import *
#from mayavi import mlab
import itertools
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex']        = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.family']       = 'sans-serif'
#mpl.rcParams['font.family']        = 'DejaVu Sans'
mpl.rcParams['lines.linewidth']    = 2
mpl.rcParams['lines.markersize']   = 6
mpl.rcParams['font.size']          = 10 
mpl.rcParams['legend.fontsize']    = 20
mpl.rcParams['axes.titlesize']     = 18
mpl.rcParams['axes.labelsize']     = 18
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
mpl.rcParams['figure.titlesize']   = 20
#mpl.rcParams['figure.figsize']     = [5,6]

import matplotlib.image as mpimg
import matplotlib.patches as patches
from matplotlib.path import Path



import os
path0 = os.getcwd()
layer, x, y, z, site = np.loadtxt('Lattice', unpack=True)    








a=0.57
Dim = 1588






######################################
######################################
######################################
######################################
cte = 1
for U in [0.5,0.7,0.9,1.1]:
        print(U)
        os.chdir('job/U'+str(cte)+'')
        site, Nup, Ndn = np.loadtxt('Dens_Nx15_Ny15_'+str(U)+'_FM.dat', unpack=True)
        Mz=(Nup-Ndn)/2

        U,EG = np.loadtxt('E_GState.dat', unpack=True)

        ######################################
        ######################################
        ######################################
        ######################################
        fig = plt.figure(figsize=(7,4))
        ax = fig.add_axes([0.02, 0.02, 0.99, 0.99])
        #ax.set_xlim(-max(x), max(x))
        #ax.set_ylim(-max(y), max(y))
        #ax.set_xticks([])
        #ax.set_yticks([])
        
        s = ax.scatter(x, y, s=50, c=Mz, alpha=0.5,
                       cmap='Reds',edgecolors='k')
        
        #s = ax.scatter(x, y, s=5000*abs(Mz), c=Mz, alpha=0.5,
        #               cmap='Reds',edgecolors='k')
        
        cbar_ax = fig.add_axes([0.85, 0.1, 0.02, 0.3])
        fig.colorbar(s, cax=cbar_ax, orientation="vertical")
        
        '''
        for i,j in itertools.product(range(Dim), range(Dim)):
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            radius2 = np.sqrt(dx*dx+dy*dy)
                
            if (z[i] == 0 and z[j] ==0) and (a-0.1 < radius2 < a+0.1):
                verts = [(x[i],y[i]), (x[j],y[j])]
                codes = [Path.MOVETO,Path.LINETO]
                path = Path(verts, codes)                                          
                ax.add_patch(patches.PathPatch(path, color='k', lw=0.05))
                
            if (z[i] > 0 and z[j] > 0) and (a-0.1 < radius2 < a+0.1):
                verts = [(x[i],y[i]), (x[j],y[j])]
                codes = [Path.MOVETO,Path.LINETO]
                path = Path(verts, codes)                                          
                ax.add_patch(patches.PathPatch(path, color='k', lw=0.05))
        '''
        
        ax.axes.text(0.86, 0.4, r'$m_z(\vec{r})$',
                     transform=ax.transAxes, 
                     rotation=0,
                     fontsize=16,
                     color='k')
        
        
        ax.axes.text(0.05, 0.95, r'$U/t='+str(U)+'$',
                     transform=ax.transAxes, 
                     rotation=0,
                     fontsize=16,
                     color='k')
        
        ax.axes.text(0.05, 0.85, r'$E_G='+str(EG)+'$',
                     transform=ax.transAxes, 
                     rotation=0,
                     fontsize=16,
                     color='k')
        
        
        ax.axis('off')
        
        cte +=1
        os.chdir(path0)



        #plt.savefig('fig_'+str(U)+'.png',dpi=300)
        plt.savefig('fig_'+str(U)+'.pdf')
        plt.show()



