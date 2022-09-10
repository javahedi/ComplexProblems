#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 11:42:29 2020

@author: javadvahedi
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
import matplotlib as mpl
mpl.rcParams['text.usetex']        = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.family']       = 'sans'
#mpl.rcParams['font.family']        = 'DejaVu Sans'
mpl.rcParams['lines.linewidth']    = 2
mpl.rcParams['lines.markersize']   = 6
mpl.rcParams['font.size']          = 10
mpl.rcParams['legend.fontsize']    = 10
mpl.rcParams['axes.titlesize']     = 10
mpl.rcParams['axes.labelsize']     = 10
mpl.rcParams['xtick.major.size']   = 3.5
mpl.rcParams['ytick.major.size']   = 3.5
mpl.rcParams['xtick.major.width']  = 1
mpl.rcParams['ytick.major.width']  = 1
mpl.rcParams['xtick.direction']    = 'in'
mpl.rcParams['ytick.direction']    = 'in'
mpl.rcParams['xtick.bottom']       = 'on'
mpl.rcParams['xtick.top']          = 'on'
mpl.rcParams['ytick.left']         = 'on'
mpl.rcParams['ytick.right']        = 'on'
mpl.rcParams['figure.titlesize']   = 10
#mpl.rcParams['figure.figsize']     = [5,6]



markers = ['o','s','v','*']
colors = ['blue','orange','green','red']
          


u0, Mmax0, tot_m0, tot_mp0 = np.loadtxt('magnetization_22-23', unpack=True)
u1, Mmax1, tot_m1, tot_mp1 = np.loadtxt('magnetization_25-26', unpack=True)
u2, Mmax2, tot_m2, tot_mp2 = np.loadtxt('magnetization_30-31', unpack=True)


plt.figure(figsize=(5,3.5))

plt.minorticks_on()

size=8
plt.plot(u2/2.7, tot_mp2 ,'o',
            label=r'(30,31), $\theta=1.08^\circ$',
            color='red',
            alpha=0.5,
            markersize=size,
            markeredgecolor='k')
            
            
plt.plot(u1/2.7, tot_mp1 ,'^',
            label=r'(25,26), $\theta=1.30^\circ$',
            color='g',
            alpha=0.5,
            markersize=size-1,
            markeredgecolor='k')
     

plt.plot(u0/2.7, tot_mp0 ,'s',
       label=r'(22,23), $\theta=1.47^\circ$',
       color='b',
       alpha=0.5,
       markersize=size-1,
       markeredgecolor='k')

"""
textstr = '\n'.join((
            r'$\theta=%.2f^\circ,~N=%i$' % (1.08,1116),
            r'moir\'e size $D=%.2f~nm$' % (12.97)
                      ))
plt.text(1., 0.3, textstr,
             fontsize=11,
             verticalalignment='top',
             color='r')
             

textstr = '\n'.join((
            r'$\theta=%.2f^\circ,~N=%i$' % (1.30,7804),
            r'moir\'e size :',
            r'$D=%.2f~nm$' % (10.85)
                      ))
plt.text(0.3, 1.7, textstr,
             fontsize=11,
             verticalalignment='top',
             color='g')
"""

plt.legend(ncol=1,loc='upper left',frameon=False,fontsize=12)

plt.xlim([0.5/2.7,5.0/2.7])
plt.ylim([-0.02,3.5])
plt.xlabel(r'$U/t_0$', fontsize = 16)
plt.ylabel(r'$M_{\rm total}$', fontsize = 16)
plt.savefig('fig2n.pdf', bbox_inches='tight')
plt.show() #; plt.close(fig)



