import numpy as np
#from mayavi.mlab import *
#from mayavi import mlab
#import itertools
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True
plt.rcParams["font.family"] = "Times New Roman"
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar

#import matplotlib.patches as patches
#from matplotlib.path import Path


a=0.57
Dim = 868
U = 1.0

layer, x, y, z, site = np.loadtxt('Lattice_868.dat', unpack=True)    
site, Nup, Ndn = np.loadtxt('Dens_Nx6_Ny6_'+str(U)+'_AF.dat', unpack=True)


Mz=(Nup-Ndn)/2

L1 = np.array([14.722431864335455,      0.50000000000000000      ])
L2 = np.array([6.9282032302755097,       13.000000000000000      ])

xnew=[]
ynew=[]
znew=[]
Mznew=[]
for n in range(-3,3):
    for m in range(-3,3):
        for i,j in zip(range(Dim),range(Dim)):
            xx, yy, zz = (m*L1+n*L2)[0], (m*L1+n*L2)[1], z[i]
            xnew.append(x[i]+xx)
            ynew.append(y[i]+yy)
            znew.append(z[i]+zz)
            Mznew.append(Mz[i])
                





fig = plt.figure(figsize=(4,4))
ax1 = fig.add_axes([0.02, 0.02, 0.96, 0.85])
ax1.set_xlim(-20, 20)
ax1.set_ylim(-20, 20)
ax1.set_xticks([])
ax1.set_yticks([])

s = ax1.scatter(xnew, ynew, s=5*np.ones(len(Mznew)), c=Mznew, alpha=1,
               cmap='bwr')#,edgecolors='k')


cbar_ax = fig.add_axes([0.2, 0.92, 0.75, 0.04])
cbar_ax.tick_params(direction='in',axis='both', color='k',
               left='on', top='on', right='on', bottom='on')
cbar_ax.xaxis.set_ticks_position("top")
fig.colorbar(s, cax=cbar_ax, orientation="horizontal")

'''
ax1_divider = make_axes_locatable(ax1)
# add an axes above the main axes.
cax1 = ax1_divider.append_axes("top", size="5%", pad="1%",axes_class="5%")
cb = colorbar(s, cax=cax1, orientation="horizontal")
# change tick position to top. Tick position defaults to bottom and overlaps
# the image.
cax1.xaxis.set_ticks_position("top")
#colorbar_ax =  fig.add_axes([0.825, 0.65, 0.03, 0.2])
#cb = fig.colorbar(s, cax=colorbar_ax)
'''
ax1.axes.text(0.02, 1.02, r'$\langle s^z_i\rangle$',
             transform=ax1.transAxes, 
             rotation=0,
             fontsize=20,
             color='k')


ax1.axes.text(0.02, 0.9, r'$U/t=1.3$',
             transform=ax1.transAxes, 
             rotation=0,
             fontsize=18,
             color='k')






ax2 = fig.add_axes([0.6, 0.03, 0.37, 0.32])
ax2.set_xticks([])
ax2.set_yticks([])
s = ax2.scatter(xnew, ynew, s=50*np.ones(len(Mznew)), c=Mznew, alpha=0.5,
               cmap='bwr',edgecolors='k')
ax2.set_xlim(-2, 2)
ax2.set_ylim(-2, 2)

#ax.axis('off')
#ax.grid('on')
#plt.tight_layout()
#plt.savefig('Mz_U_'+str(U)+'_new2.png',dpi=200)#,transparent=True)
#plt.savefig('Mz_U_'+str(U)+'_new.pdf')

#plt.savefig('Mz_U_'+str(U)+'_MFT.png',dpi=200)
plt.show()
