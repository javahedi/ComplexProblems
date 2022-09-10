
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path



def TwoDimPlot(Dim,a):
    site, x, y, z, ABCD = np.loadtxt('Lattice_'+str(Dim)+'.dat', unpack=True)

    N = len(site)
    Max_x = np.ceil(max(x))
    Max_y = np.ceil(max(y))
    Min_x = np.floor(min(x))
    Min_y = np.floor(min(y))

    
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, aspect='equal')
    
    ax.axis([Min_x-a, Max_x+a, Min_y-a, Max_y+a])
    for i, j, s in zip(range(N), range(N), site):
        if z[i]>0:
            circle = patches.Circle((x[i], y[j]), radius = 0.1*a , color = 'blue' )
        else:
            circle = patches.Circle((x[i], y[j]), radius = 0.1*a , color = 'red' )

            
        ax.add_patch(circle)
    
        ax.text(x[i], y[j], str(int(s)), color="black", fontsize=15)
        
    rectangle = patches.Rectangle((Min_x-2, -23),         # (x,y)
                                      Max_x+a+2,          # width
                                      33,          # height
                                      #fill=False,
                                      facecolor="#00ffff",
                                      alpha = 0.5,
                                      linestyle='dashed')
    ax.add_patch(rectangle)
    #http://matthiaseisen.com/pp/patterns/p0203/


    # lattice parameters
    r1 = 22.40  # nm
    r2 = 22.80
    r1x, r1y, r1z = 15.03, 16.60,   0.0
    r2x, r2y, r2z =  7.86,   0.0, 21.40
    a = 45.80
    b = 33.20
    
    r3 = np.sqrt((r1x+2*r2x)**2+r1y**2)
    r4 = np.sqrt((r1x+r2x)**2+r1y**2)
    r5 = 2*r1x+r2x
    
    for i in range(N):
        for j in range(N):        
            dx = x[i] - x[j]
            dy = y[i] - y[j] 
            dz = z[i] - z[j]
                
            if abs(dz) == 0 and r1-0.1 < np.sqrt(dx*dx+dy*dy) < r1+0.1: # t1 up or bottom 
                if z[i]> 0:
                    verts = [(x[i],y[i]), (x[j],y[j])]
                    codes = [Path.MOVETO,Path.LINETO]
                    path = Path(verts, codes)
                    plt.gca().add_patch(patches.PathPatch(path, color='red', lw=0.5))
                else:
                    verts = [(x[i],y[i]), (x[j],y[j])]
                    codes = [Path.MOVETO,Path.LINETO]
                    path = Path(verts, codes)
                    plt.gca().add_patch(patches.PathPatch(path, color='blue', lw=0.5))
            if abs(dz) > 0 and r2-0.1 < np.sqrt(dx*dx+dy*dy+dz*dz) < r2+0.1: # t2 up <-> bottom
                verts = [(x[i],y[i]), (x[j],y[j])]
                codes = [Path.MOVETO,Path.LINETO]
                path = Path(verts, codes)
                plt.gca().add_patch(patches.PathPatch(path, color='black', lw=1.5))
            '''    
            if abs(dz) == 0 and r3-1 < np.sqrt(dx*dx+dy*dy) < r3+1: # t3 up or bottom
                if color[i] == 0:
                    verts = [(x[i],y[i]), (x[j],y[j])]
                    codes = [Path.MOVETO,Path.LINETO]
                    path = Path(verts, codes)
                    plt.gca().add_patch(patches.PathPatch(path, color='red', lw=0.5))
                else:
                    verts = [(x[i],y[i]), (x[j],y[j])]
                    codes = [Path.MOVETO,Path.LINETO]
                    path = Path(verts, codes)
                    plt.gca().add_patch(patches.PathPatch(path, color='blue', lw=0.5))
                  
            if abs(dz) > 0 and r4-1 < np.sqrt(dx*dx+dy*dy) < r4+1: # t2 up <-> bottom
                verts = [(x[i],y[i]), (x[j],y[j])]
                codes = [Path.MOVETO,Path.LINETO]
                path = Path(verts, codes)
                plt.gca().add_patch(patches.PathPatch(path, color='black', lw=1.5))
                
            if abs(dz) > 0 and r5-1 < np.sqrt(dx*dx+dy*dy) < r5+1: # t2 up <-> bottom
                verts = [(x[i],y[i]), (x[j],y[j])]
                codes = [Path.MOVETO,Path.LINETO]
                path = Path(verts, codes)
                plt.gca().add_patch(patches.PathPatch(path, color='black', lw=1.5))
            '''        
#    print('r1',r1)
#    print('r2',r2)
#    print('r3',r3)
#    print('r4',r4)
#    print('r5',r5)

    ax.axis('off')
    plt.savefig('Lattice.pdf')   
      
    return plt.show()

