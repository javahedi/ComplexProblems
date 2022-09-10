import numpy as np
#from mayavi.mlab import *
#from mayavi import mlab
import itertools
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True
#from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.patches as patches
from matplotlib.path import Path



def Spatial(Nx,Ny,Size,Dim,U):
    
    site, Nup, Ndn = np.loadtxt('Density_U_'+str(U)+'.txt', unpack=True)
    site, x, y, z, label = np.loadtxt('Lattice_'+str(Size)+'.txt', unpack=True)    
    
    
    
    Mz=np.subtract(Nup,Ndn)/2
    
    
    fig = plt.figure(figsize=(4,6))
    ax = fig.add_axes([0.01, 0.01, 0.95, 0.95])
    #plt.axis([-max(x), max(x), -max(y), max(y)])
    ax.set_xlim(-30, max(x)+35)
    ax.set_ylim(-10, y[Dim])
    
    
    s=ax.scatter(x[:Dim], y[:Dim], s=1000*abs(Mz), c=Mz, alpha=1.0,
               cmap='bwr',edgecolors='k')

    colorbar_ax =  fig.add_axes([0.86, 0.2, 0.02, 0.6])
    fig.colorbar(s, cax=colorbar_ax)
    
    
    
    # lattice parameters
    r1 = 22.40  # nm
    r2 = 22.80
    r1x, r1y, r1z = 15.03, 16.60,   0.0
    r2x, r2y, r2z =  7.86,   0.0, 21.40
    
    
    
    for i in range(Dim):
        for j in range(Dim):        
            dx = x[i] - x[j]
            dy = y[i] - y[j] 
            dz = z[i] - z[j]
                
            if abs(dz) == 0 and r1-0.1 < np.sqrt(dx*dx+dy*dy) < r1+0.1: # t1 up or bottom 
                if z[i] > 0:
                    verts = [(x[i],y[i]), (x[j],y[j])]
                    codes = [Path.MOVETO,Path.LINETO]
                    path = Path(verts, codes)
                    ax.add_patch(patches.PathPatch(path, color='k', lw=0.1))
                else:
                    verts = [(x[i],y[i]), (x[j],y[j])]
                    codes = [Path.MOVETO,Path.LINETO]
                    path = Path(verts, codes)
                    ax.add_patch(patches.PathPatch(path, color='k', lw=0.1))
            if abs(dz) > 0 and r2-0.1 < np.sqrt(dx*dx+dy*dy+dz*dz) < r2+0.1: # t2 up <-> bottom
                verts = [(x[i],y[i]), (x[j],y[j])]
                codes = [Path.MOVETO,Path.LINETO]
                path = Path(verts, codes)
                ax.add_patch(patches.PathPatch(path, color='k', lw=0.1))
    
    

    textstr = '\n'.join((    
        r'$U/|t_1|=%.2f$' % (U),
        r'$Nx,Ny=%i,%i$' % (Nx,Ny//2),
                              ))
                        
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    
    ax.axes.text(0.35, 0.99, textstr, transform=ax.transAxes, fontsize=14,
           verticalalignment='top', bbox=props)
    
    
    ax.axis('off')
    
    
    #plt.savefig('Mz_Nx_'+str(Nx)+'_Ny_'+str(Ny//2)+'_U_'+str(U)+'_.pdf')
    plt.savefig('Mz_Nx_'+str(Nx)+'_Ny_'+str(Ny//2)+'_U_'+str(U)+'_.png')


    return #plt.show()
