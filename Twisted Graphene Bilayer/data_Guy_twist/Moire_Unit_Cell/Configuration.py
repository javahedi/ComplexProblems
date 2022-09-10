#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 14:40:33 2018

@author: javad
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import itertools

import matplotlib as mpl
mpl.rcParams['text.usetex']        = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.family']        = 'sans'
#mpl.rcParams['font.family']         = 'sans'
mpl.rcParams['lines.linewidth']    = 2
mpl.rcParams['lines.markersize']   = 6
mpl.rcParams['font.size']          = 18 
mpl.rcParams['legend.fontsize']    = 20
mpl.rcParams['axes.titlesize']     = 14
mpl.rcParams['axes.labelsize']     = 18
mpl.rcParams['xtick.major.size']   = 0
mpl.rcParams['ytick.major.size']   = 0
mpl.rcParams['xtick.major.width']  = 0
mpl.rcParams['ytick.major.width']  = 0
mpl.rcParams['xtick.direction']    = 'in'
mpl.rcParams['ytick.direction']    = 'in'
mpl.rcParams['xtick.bottom']       = True
mpl.rcParams['xtick.top']          = True
mpl.rcParams['ytick.left']         = True
mpl.rcParams['ytick.right']        = True
mpl.rcParams['figure.titlesize']   = 20
#mpl.rcParams['figure.figsize']     = [6,5]


def TwoDimPlot(r,m0,a,d,R1,R2,theta,N,xmax,ymax):
    site, x, y, z=  np.loadtxt('Lattice_r_'+str(r)+'m0_'+str(m0)+'.dat', unpack=True)

    n = len(site)
    
    fig = plt.figure(figsize=(5,4))
    ax = fig.add_axes([0.05, 0.01, 0.9, 0.9],facecolor='w')
        
        
    #uncomment if want to add rhomboid unit ceel at the center of Figure
    
    AB = 1*(R1+R2)/2
    
    
    verts = [(0,0)-AB, (R1[0],R1[-1])-AB]
    codes = [Path.MOVETO,Path.LINETO]
    path = Path(verts, codes)                                          
    plt.gca().add_patch(patches.PathPatch(path, color='black', lw=1)) 
    
    verts = [(0,0)-AB, (R2[0],R2[-1])-AB]
    codes = [Path.MOVETO,Path.LINETO]
    path = Path(verts, codes)                                          
    plt.gca().add_patch(patches.PathPatch(path, color='black', lw=1)) 
    
    
    verts = [(R1[0],R1[-1])-AB, (R1[0]+R2[0],R1[-1]+R2[-1])-AB]
    codes = [Path.MOVETO,Path.LINETO]
    path = Path(verts, codes)                                          
    plt.gca().add_patch(patches.PathPatch(path, color='black', lw=1)) 
    
    
    verts = [(R2[0],R2[-1])-AB, (R1[0]+R2[0],R1[-1]+R2[-1])-AB]
    codes = [Path.MOVETO,Path.LINETO]
    path = Path(verts, codes)                                          
    plt.gca().add_patch(patches.PathPatch(path, color='black', lw=1)) 
    
    
    
    #uncomment if connection between sites are interested
    '''
    for i,j in itertools.product(range(n), range(n)):
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            radius2 = np.sqrt(dx*dx+dy*dy)
                
            if (z[i] == 0 and z[j] ==0) and (0.55 < radius2 < 0.58):
                verts = [(x[i],y[i]), (x[j],y[j])]
                codes = [Path.MOVETO,Path.LINETO]
                path = Path(verts, codes)                                          
                plt.gca().add_patch(patches.PathPatch(path, color='b', lw=0.5))  
                      
            if (z[i] > 0 and z[j] > 0) and (0.55 < radius2 < 0.58):
                verts = [(x[i],y[i]), (x[j],y[j])]
                codes = [Path.MOVETO,Path.LINETO]
                path = Path(verts, codes)                                          
                plt.gca().add_patch(patches.PathPatch(path, color='r', lw=0.5))               
  
    '''
    
    
    
    x_b,x_t = x[:n//2],x[n//2:]
    y_b,y_t = y[:n//2],y[n//2:]
    ax.scatter(x_b,y_b,s=4,color='r')
    ax.scatter(x_t,y_t,s=4,color='b')

   

    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_ticklabels([])
    
    ax.set_ylim([-ymax, ymax])
    ax.set_xlim([-xmax, xmax])
    
    ax.set_title(r'(a)',loc='left')
    plt.subplots_adjust(left=0.01, bottom=0.1, right=0.99, top=0.9, wspace=0.0, hspace=0.0)
    plt.savefig('fig1a.pdf',dpi=300, bbox_inches='tight')
    plt.savefig('fig1a.png',dpi=600, bbox_inches='tight')
    return plt.show()

