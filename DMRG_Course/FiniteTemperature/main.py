import numpy as np
from copy import deepcopy
import tensornetwork as tn
import tqdm
import time
import os, sys
os.environ['MKL_NUM_THREADS']='4' # set number of MKL threads to run in parallel
#os.environ['OMP_NUM_THREADS']='1' # set number of OpenMP threads to run in parallel


###############################################
def CV():
    ######## Define the simulation parameter ######################
    L      = 20
    d      = 4
    chiMax = 40
    BC     = 'OBC'
    ######## Define the model parameter ###########################
    Jxy    = 1.0
    Jz     = 0.0
    hz     = 0.
    hx     = 0.0
    delta  = 0.1
    tfin   = 8.0
    N = int(tfin/delta)
    ######### initial MPS state ##########
    #l, G, chi  = tn.initialMat(2*L, 2, chiMax, BC)
    

    #print(G[0].real)
    #print("=="*10)
    #print(G[1].real)
    #print("=="*10)
    #print(np.squeeze(np.tensordot(G[0],G[1],axes=(2,1)).real))
    
    ######### Mix  MPS state ##########
    #AA, BB = tn.lG_AABB(l, G)
    ######### Back to  vidal ##########
    #l, G, norm = tn.MPS_Vidal(AA)
    
    #print("=="*10)
    #print(np.squeeze(np.tensordot(BB[0],BB[1],axes=(2,1)).real))
    #print("=="*10)
    
    l, G, chi  = tn.initialMat_mix(L, d, chiMax, BC)
    #print(np.squeeze(np.tensordot(G[0],G[1],axes=(2,1)).real))
   
    ######### normalize MPS  ##########
    l, G, norm = tn.normalize_left (l,  G)
    l, G, norm = tn.normalize_right(l,  G)
    
   
    
    H = tn.hamiltinian(Jxy,Jz,hz,hx,True)

    
    
    f1 = open('cv.txt','w')
    f2 = open('energy.txt','w')
    ee = np.zeros(N)
    tem = np.zeros(N)
    for i in tqdm.tqdm(range(N)):
    #for i in range(N):
        l, G, disc = tn.time_evolution(l, G, H, BC, -delta/2, 1)
        normalize  = np.squeeze(tn.overlap(l, G, l, G).real).item()
        l0 = deepcopy(l)
        G0 = deepcopy(G)
        
        E = 0.
        for pos in range(L-1):
            E += tn.expvalue2(l,G,H ,pos).real
        E /=normalize
        ee[i]=E
            
        l0, G0, disc = tn.evolve_bonds(l0,G0,H,np.arange(0,L))
        E2  = tn.overlap(l0, G0, l0, G0).real
        E2 /=normalize
        
        T = (i+1)*delta
        tem[i]=T
        cv = np.squeeze((E2-E**2)/T**2).item()
        f1.write('%f    %f\n'%(T,cv/L))
        f2.write('%f    %f   %f\n'%(T,E/L,E2/L))
        f1.flush()
        f2.flush()
    f1.close()
    f2.close()
    
    cvv = np.diff(ee)/tem[:-1]
    np.savetxt('cvv2.txt',cvv/L)
    
CV()
