import numpy as np
from copy import deepcopy
import tensornetwork as tn
import tqdm
import time
import os, sys
os.environ['MKL_NUM_THREADS']='2' # set number of MKL threads to run in parallel
#os.environ['OMP_NUM_THREADS']='1' # set number of OpenMP threads to run in parallel




###############################################
def ground_state():
    ######## Define the simulation parameter ######################
    L      = 100
    d      = 2
    chiMax = 100
    BC     = 'OBC'
    ######## Define the model parameter ###########################
    Jxy    = 0.0
    Jz     = -1.0
    hz     = 0.0
    hx     = 1.5
    tfin   = 20.0
    ######### initial MPS state ##########
    l, G, chi  = tn.initialMat(L, d, chiMax, BC)
    l, G, norm = tn.normalize_left (l,  G)
    #l, G, norm = tn.normalize_right
    
    H = tn.hamiltinian(Jxy,Jz,hz,hx)
    print("Loop over temperatue to find ground state")
    for delta in [0.1]:
       N = int(tfin/delta)
       for i in tqdm.tqdm(range(N)):
           l, G, disc = tn.time_evolution(l, G, H, BC, -delta, 1)
    print(" Ground state found -:) ")
    return l, G
    


###############################################
def quench_dynamic(l,G):
    l0, G0 = deepcopy(l), deepcopy(G)
    ######## Define the simulation parameter ######################
    L      = 100
    d      = 2
    chimax = 100
    BC     = 'OBC'
    ######## Define the model parameter ###########################
    Jxy    = 0.0
    Jz     = -1.0
    hz     = 0.
    hx     = 0.2
    delta  = 0.1
    tfin   = 5.0
    N      = int(tfin/delta)
    ######### initial MPS state ##########
    H = tn.hamiltinian(Jxy,Jz,hz,hx)
    
    f1 = open('quench_PRB_fig1.txt','w')
    disc_tot = 0
    for i in tqdm.tqdm(range(N)):
          l, G, disc = tn.time_evolution(l, G, H, BC, -1J*delta, 1)
          OL = np.abs(np.squeeze(tn.overlap(l0,G0, l,G)).item())
          rt = -np.log(np.abs(OL)**2)/L
          f1.write(str((i+1)*delta)+"\t"+str(OL.real)+"\t"+str(rt.real)+"\n")
          f1.flush()
          #print((i+1)*delta, disc)
          disc_tot +=disc
    print(' total discar : ',disc_tot)

    
l, G = ground_state()
print(" ")
quench_dynamic(l,G)

