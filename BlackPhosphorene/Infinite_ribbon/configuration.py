
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path



def TwoDimPlot(Dim,a):
    site, x, y, z, ABCD = np.loadtxt('Lattice_'+str(Dim)+'.txt', unpack=True)

    N = len(site)
    Max_x = np.ceil(max(x))
    Max_y = np.ceil(max(y))
    Min_x = np.floor(min(x))
    Min_y = np.floor(min(y))

    
    #fig = plt.figure(figsize=(10,4))
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    
    ax.axis([Min_x-a, Max_x+a, Min_y-a, Max_y+a])
    for i, j, s in zip(range(N), range(N), site):
        if ABCD[i]==0 or ABCD[i]==2:
            circle = patches.Circle((x[i], y[j]), radius = 0.1*a , color = 'blue' )
        else:
            circle = patches.Circle((x[i], y[j]), radius = 0.1*a , color = 'red' )

            
        ax.add_patch(circle)
    
        ax.text(x[i], y[j], str(int(s)), color="black", fontsize=12)
        
    rectangle = patches.Rectangle((Min_x-3, -10),         # (x,y)
                                      Max_x+a+3,          # width
                                      32,          # height
                                      #fill=False,
                                      facecolor="#00ffff",
                                      alpha = 0.5,
                                      linestyle='dashed')
    ax.add_patch(rectangle)
    #http://matthiaseisen.com/pp/patterns/p0203/

    for i in range(N):
        for j in range(N):        
            d_x = x[i] - x[j]
            d_y = y[i]- y[j]
            radius2 = d_x*d_x+d_y*d_y
            if a-0.1 < radius2 < a+0.1:
                verts = [(x[i],y[i]), (x[j],y[j])]
                codes = [Path.MOVETO,Path.LINETO]
                path = Path(verts, codes)
                plt.gca().add_patch(patches.PathPatch(path, color='black', lw=0.5))

    #ax.axis('off')
    plt.savefig('Lattice.pdf')   
      
    return plt.show()

