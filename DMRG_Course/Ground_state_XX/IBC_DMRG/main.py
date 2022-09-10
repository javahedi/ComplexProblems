import numpy as np
from copy import deepcopy
import tensornetwork as tn
import tqdm
import time
import os, sys
os.environ['MKL_NUM_THREADS']='4' # set number of MKL threads to run in parallel
#os.environ['OMP_NUM_THREADS']='1' # set number of OpenMP threads to run in parallel


###############################################
def ground_state():
    ######## Define the simulation parameter ######################
    L      = 2
    d      = 2
    chiMax = 100
    BC     = 'IBC'
    ######## Define the model parameter ###########################
    Jxy    = 1.0
    Jz     = 0.0
    hz     = 0.
    hx     = 0.0
    tfin   = 40.0
    ######### initial MPS state ##########
    l, G, chi  = tn.initialMat(L, d, chiMax, BC)
    #l, G, norm = tn.normalize_left (l,  G)
    #l, G, norm = tn.normalize_right(l,  G)
    
    H = tn.hamiltinian(Jxy,Jz,hz,hx)
    for delta in [0.01]:
       N = int(tfin/delta)
       for i in tqdm.tqdm(range(N)):
           l, G, disc = tn.time_evolution(l, G, H, BC, -delta, 1)
           #print(len(tn.overlap(l, G, l, G)))
       sum = 0.
       for pos in range(L):
           sum += tn.expvalue2_inf(l,G,H,pos).real
       
       print(sum/L,"(DMRG with dt : ", delta,")")
       print("     ")
    print("(E0/L =", -1.0/np.pi)
    
    
ground_state()
