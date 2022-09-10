import numpy as np
#from mayavi.mlab import *
#from mayavi import mlab
import itertools
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True
plt.rcParams["font.family"] = "Times New Roman"

import matplotlib.patches as patches
from matplotlib.path import Path


a=1.42
Dim = 868
m,n=8,9
U = 1.0
site, x, y, z = np.loadtxt('unit_Cell_m_'+str(m)+'_n_'+str(n)+'.dat', unpack=True)    
#site, Nup, Ndn = np.loadtxt('Density_U_'+str(U)+'_.dat', unpack=True)
site, Nup, Ndn = np.loadtxt('Dens_Nx6_Ny6_'+str(U)+'_AF.dat', unpack=True)



Mz=np.subtract(Nup,Ndn)/2


fig = plt.figure()
ax = fig.add_axes([0.01, 0.01, 0.95, 0.95])
#plt.axis([-max(x), max(x), -max(y), max(y)])
ax.set_xlim(-1, max(x))
ax.set_ylim(-1, max(y))

s = ax.scatter(x, y, s=10000*abs(Mz), c=Mz, alpha=1.0,
               cmap='bwr',edgecolors='k')


colorbar_ax =  fig.add_axes([0.85, 0.05, 0.02, 0.6])
#colorbar_ax =  fig.add_axes([0.05, 0.02, 0.8, 0.04])

fig.colorbar(s, cax=colorbar_ax)
#cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=20)

for i,j in itertools.product(range(Dim), range(Dim)):
    dx = x[i] - x[j]
    dy = y[i] - y[j]
    radius2 = np.sqrt(dx*dx+dy*dy)
        
    if (z[i] == 0 and z[j] ==0) and (a-0.1 < radius2 < a+0.1):
        verts = [(x[i],y[i]), (x[j],y[j])]
        codes = [Path.MOVETO,Path.LINETO]
        path = Path(verts, codes)                                          
        ax.add_patch(patches.PathPatch(path, color='k', lw=0.1))
        
    if (z[i] > 0 and z[j] > 0) and (a-0.1 < radius2 < a+0.1):
        verts = [(x[i],y[i]), (x[j],y[j])]
        codes = [Path.MOVETO,Path.LINETO]
        path = Path(verts, codes)                                          
        ax.add_patch(patches.PathPatch(path, color='k', lw=0.1))
              
   


'''
ax.axes.text(-10.0, 4.0, r'$\theta_{eff}=1.08^\circ$', color='k', fontsize=20,
              bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
ax.axes.text(-9.0, 3.0, r'$U='+str(U)+'~eV$', color='k', fontsize=15,
              bbox={'facecolor':'red', 'alpha':0.5, 'pad':10},
              rotation=60)
'''

textstr = '\n'.join((
    r'$\theta_{eff}=%.2f^\circ$' % (1.08),    
    r'$U=%.2f~(eV)$' % (U),
    r'$Nx=Ny=%i$' % (6),
                          ))
                    
props = dict(boxstyle='round', facecolor='wheat', alpha=0.1)

ax.axes.text(-0.0, 0.82, textstr, transform=ax.transAxes, fontsize=18,
        verticalalignment='top', bbox=props)

ax.axes.text(0.02, 0.9, r'(b)',transform=ax.transAxes, fontsize=18)

#ax.set_title(r'(b)',loc='left',fontsize=18)
ax.axis('off')
#plt.tight_layout()
plt.savefig('Mz_NxNy6_'+str(U)+'.png',transparent=True,dpi=200)
plt.show()
