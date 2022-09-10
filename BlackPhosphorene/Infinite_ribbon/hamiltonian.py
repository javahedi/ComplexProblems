from math import sqrt
import numpy as np

def Hamil(t,e0,filename):

    
    site, x, y, z, ABCD = np.loadtxt(filename, unpack=True)   
    N=len(site)     # size of sytem    
    # lattice parameters
    r1 = 22.40  # nm
    r2 = 22.80
    r1x, r1y, r1z = 15.03, 16.60,   0.0
    r2x, r2y, r2z =  7.86,   0.0, 21.40    
    #a = 45.80
    #b = 33.20
    
    r3 = sqrt((r1x+2*r2x)**2+r1y**2)
    r4 = sqrt((r1x+r2x)**2+r1y**2)
    r5 = 2*r1x+r2x    
    Ham=np.zeros((N,N))    
   
    for i in range(N):
        if z[i] > 0: 
            Ham[i][i] = e0/2
        else:
            Ham[i][i] = -e0/2          
            
        #cte = 0
        for k in range(N):
            dx = x[k] - x[i]
            dy = y[k] - y[i]
            dz = z[k] - z[i]
            if abs(dz) == 0 and r1-0.1 < sqrt(dx*dx+dy*dy) < r1+0.1: # t1 up or bottom
                Ham[int(site[i])-1][int(site[k])-1] = t[0]
            if abs(dz) > 0 and  r2-0.1 < sqrt(dx*dx+dy*dy+dz*dz) < r2+0.1: # t2 up <-> bottom
                Ham[int(site[i])-1][int(site[k])-1] = t[1]                
            if abs(dz) == 0 and r3-1 < sqrt(dx*dx+dy*dy) < r3+1: # t3 up or bottom                
                Ham[int(site[i])-1][int(site[k])-1] = t[2]               
            if abs(dz) > 0 and  r4-1 < sqrt(dx*dx+dy*dy) < r4+1: # t4 up <-> bottom               
                Ham[int(site[i])-1][int(site[k])-1] = t[3]               
            if abs(dz) > 0 and  r5-1 < sqrt(dx*dx+dy*dy) < r5+1: # t5 up <-> bottom                
                Ham[int(site[i])-1][int(site[k])-1] = t[4]
    
    
    return Ham
            
    
