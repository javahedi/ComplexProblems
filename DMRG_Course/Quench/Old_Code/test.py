import numpy as np
from copy import deepcopy
import tensornetwork as tn
from ed_ising import ising_gs
import tqdm
import time
import os, sys
os.environ['MKL_NUM_THREADS']='4' # set number of MKL threads to run in parallel
#os.environ['OMP_NUM_THREADS']='1' # set number of OpenMP threads to run in parallel


###############################################
def clean_quench():
    ######## Define the simulation parameter ######################
    L      = 20
    d      = 2
    chiMax = 100
    BC     = 'obc'
    op     = np.array([[1.,0.],[0.,-1.]])/2.
    ######## Define the model parameter ###########################
    Jxy    = 1.0
    Jz     = 0.5
    hz     = 0.
    hx     = 0.
    delta  = 0.1
    tfin   = 10.0
    N      = int(tfin/delta)
    ######### initial MPS state ##########
    l, G, chi  = tn.initialMat(L, d, chiMax, BC)
    l, G, norm = tn.normalize_left (l,  G)
    l, G, norm = tn.normalize_right(l,  G)    

    H = tn.hamiltinian(Jxy,Jz,hz,hx)

    f1 = open('mz_L20_i10.txt','w')
    f2 = open('EE_L20_i10.txt','w')
    f3 = open('error_L20_i10.txt','w')
    for i in tqdm.tqdm(range(N)):
          l, G, disc = tn.time_evolution(l, G, H, BC, -1J*delta, 1)
          
          sz  =  tn.expvalue1(l,G,op,L//2)
          ee  =  tn.entanglement_entropy(l)
          
          f1.write(str((i+1)*delta)+"\t"+str(2.0*sz.real)+"\n")
          f2.write(str((i+1)*delta)+"\t"+str(ee[L//2])+"\n")
          f3.write(str(disc)+"\n")
          f1.flush()
          f2.flush()
          f3.flush()
          #print((i+1)*delta, disc)
    f1.close()
    f2.close()
    f3.close()

clean_quench()
