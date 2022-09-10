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
          


uu, m1 ,m2, m3 = np.loadtxt('magnetization_Mono_t0.75eV_k27x27', unpack=True)

U3, Mz3 ,Mz_ave3= np.loadtxt('Mz_75_zero.dat', unpack=True)
U8, Mz8 ,Mz_ave4= np.loadtxt('Mz_75_AA.dat', unpack=True)
U9, Mz9 ,Mz_ave5= np.loadtxt('Mz_75_AB.dat', unpack=True)


u1, Mmax1, tot_m1, tot_mp1 = np.loadtxt('magnetization_ThetaEff_1.08', unpack=True)
u2, Mmax2, tot_m2, tot_mp2 = np.loadtxt('magnetization_ThetaEff_1.30', unpack=True)
u3, Mmax3, tot_m3, tot_mp3 = np.loadtxt('magnetization_ThetaEff_1.47', unpack=True)



fig, subfig = plt.subplots(2,1,figsize=(5,6))

subfig[0].minorticks_on()
subfig[1].minorticks_on()


rect = Rectangle((2.1, 0.0),   # (x,y)
                  0.2,          # width
                  0.2,          # height
                  fill=True,
                  zorder=3,
                  color='lightgrey',
                  alpha=0.6)
#subfig[1].add_patch(rect)



size=3

subfig[0].plot(u1/0.75, tot_mp1 ,
            label=r'$\theta_{\rm eff}=1.08^\circ$',
            color='red',
            marker=markers[0],
            markersize=size,
            markeredgecolor=None)
subfig[0].plot(u2/0.9, tot_mp2 ,
            label=r'$\theta_{\rm eff}=1.30^\circ$',
            color='g',
            marker=markers[1],
            markersize=size,
            markeredgecolor=None)
subfig[0].plot(u3/1.02, tot_mp3 ,
            label=r'$\theta_{\rm eff}=1.47^\circ$',
            color='blue',
            marker=markers[2],
            markersize=size,
            markeredgecolor=None)
          




subfig[0].set_xlim([0,2.5])
subfig[0].set_ylim([0,100])
subfig[0].set_ylabel(r'$M_{\rm total}$', fontsize = 18)
subfig[0].text(-0.16, 1, r'(a)',
             transform=subfig[0].transAxes,
             fontsize=16,
             color='k',
             rotation=0)




#ax0 = fig.add_axes([0.12, 0.58, 0.5, 0.28],facecolor='w')
ax0 = subfig[0].inset_axes([0.02, 0.1, 0.6, 0.4])
ax0.yaxis.tick_right()
#ax0.xaxis.tick_top()
ax0.xaxis.set_ticks_position('top') # the rest is the same
ax0.tick_params(labelbottom=False,labeltop=True)

#ax0.set_xlabel(r'$U/t^\prime$', fontsize = 14)

ax0.plot(u1/0.75, tot_mp1 ,'o',
            label=r'$\theta_{eff}=1.08^\circ$',
            color='red',
            markersize=size,
            markeredgecolor=None)
ax0.plot(u2/0.9, tot_mp2,'s',
            label=r'$\theta_{eff}=1.30^\circ$',
            color='g',
            markersize=size,
            markeredgecolor=None)
ax0.plot(u3/1.02, tot_mp3 ,'^',
            label=r'$\theta_{eff}=1.47^\circ$',
            color='blue',
            markersize=size,
            markeredgecolor=None)
          
ax0.set_xlim([0.1,1.5])
ax0.set_ylim([0,2.0])


ax0.minorticks_on()


#subfig[1].plot(U3/0.75, Mz3,
#                '-',lw=3,
#                label=r'Graphene',
#                color='orange')
subfig[1].plot(uu/0.75, m1,
                '-',lw=3,
                label=r'Graphene',
                 color='orange')



subfig[1].plot(u1/0.75, Mmax1,'o',
            color='red',
            markersize=size,
            markeredgecolor=None)
subfig[1].plot(u2/0.9, Mmax2 ,'s',
            color='g',
            markersize=size,
            markeredgecolor=None)
subfig[1].plot(u3/1.02, Mmax3, '^',
            color='blue',
            markersize=size,
            markeredgecolor=None)


subfig[1].set_yticks([0.0,0.05,0.1,0.15,0.2])
subfig[1].set_yticklabels([0.0,0.05,0.1,0.15,0.2])



x00,x01,x02,y00,y01,y02,=0.25,0.32,1.0,0.07,0.07,0.09
x10,x11,x12,y10,y11,y12=0.25,0.32,1.0,0.0 ,0.0 ,0.0
subfig[1].quiver([x00,x01,x02],
                 [y00,y01,y02],
                 [x10-x00,  x11-x01, x12-x02],
                 [y10-y00,  y11-y01, y12-y02],
                angles='xy', scale_units='xy', scale=1,
                zorder=3, color=['g','r','b'], alpha=0.5,
                width=0.01, headwidth=3., headlength=4.)
       
subfig[1].text(0.05, 0.15,
                r'$U_{c1}/t_0^\prime\approx 0.21$',
                transform=subfig[1].transAxes,
                fontsize=8,
                color='g',
                rotation=90)


subfig[1].text(0.15, 0.15,
                r'$U_{c1}/t_0^\prime\approx 0.32$',
                transform=subfig[1].transAxes,
                fontsize=8,
                color='r',
                rotation=90)

subfig[1].text(0.42, 0.22,
                r'$U_{c1}/t_0^\prime\approx 1.0$',
                transform=subfig[1].transAxes,
                fontsize=8,
                color='b',
                rotation=90)


#subfig[1].text(0.65, 0.6, r'$U_{c2}/t^\prime\approx 2.1$',
#             transform=subfig[1].transAxes, fontsize=8,
#             color='b', rotation=90)


subfig[1].legend(ncol=1,loc='upper left',frameon=False,fontsize=12)
subfig[1].set_xlabel(r'$U/t_0^\prime$', fontsize = 18)
subfig[1].set_ylabel(r'$M_{\rm max}$', fontsize = 18)
subfig[1].set_xlim([0,2.5])
subfig[1].set_ylim([0,0.2])

subfig[1].text(-0.16, 1, r'(b)',
                transform=subfig[1].transAxes,
                fontsize=16,
                color='k', rotation=0)
                

subfig[0].legend(ncol=1,
                loc='upper left',
                frameon=False,fontsize=12)

plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.99,
                    top=0.95,
                    wspace=0.0,
                    hspace=0.12)
                    
                    


#subfig[1].text(0.85, 0.6,
#            r'$U_{c2}/t_0^\prime$',
#            transform=subfig[1].transAxes,
#            fontsize=14,
#            color='k',
#            rotation=90)
                    
plt.savefig('fig2.pdf', bbox_inches='tight')
plt.show(fig) #; plt.close(fig)



