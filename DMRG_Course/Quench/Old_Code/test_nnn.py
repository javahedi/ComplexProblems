import numpy as np
from copy import deepcopy
import tensornetwork as tn
from ed_ising import ising_gs
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
    chiMax = 50
    BC     = 'obc'
    ######## Define the model parameter ###########################
    Jxy    = 0.0
    Jz     = -1.0
    hz     = 0.
    hx     = 1.3
    delta  = 0.01
    tfin   = 5.0
    N      = int(tfin/delta)
    ######### initial MPS state ##########
    l, G, chi  = tn.initialMat(L, d, chiMax, BC)
    l, G, norm = tn.normalize_left (l,  G)
    l, G, norm = tn.normalize_right(l,  G)

    H1 = tn.hamiltinian(Jxy,Jz,hz,hx)
    H2 = H1

    disc_tot = 0
    for i in tqdm.tqdm(range(N)):
          l, G, disc = tn.time_evolution(l, G, H1, BC, -delta, 1)  # expH(i,j+2)
          
          # for next nearest neighbeour sites
          # swaping  gape  expH(i,j+2) = S(i+1,j+2) expH(i,j+1) S(i+1,j+2)
          G          = tn.swape_gate(G)                            # S(i+1,j+2)
          l, G, disc = tn.time_evolution(l, G, H1, BC, -delta, 1)  # expH(i,j+2)
          G          = tn.swape_gate(G)                            # S(i+1,j+2)
          
          
          disc_tot +=disc
    print(' total discar : ',disc_tot)
    return l, G
      
          



###############################################
def quench(l,G):
    l0, G0 = deepcopy(l), deepcopy(G)

    ######## Define the simulation parameter ######################
    L      = 100
    d      = 2
    chiMax = 50
    BC     = 'obc'
    ######## Define the model parameter ###########################
    Jxy    = 0.0
    Jz     = -1.0
    hz     = 0.0
    hx     = 1.3
    delta  = 0.1
    tfin   = 5.0
    N      = int(tfin/delta)
    
    
    H1 = tn.hamiltinian(Jxy,Jz,hz,hx)
    H2 = tn.hamiltinian(Jxy,0.6*Jz,hz,hx)

    f1 = open('quench_PRB_L_'+str(L)+'_new.txt','w')
    f2 = open('error_PRB_L_'+str(L)+'_new.txt','w')
    disc_tot = 0
    for i in tqdm.tqdm(range(N)):
          l, G, disc = tn.time_evolution(l, G, H, BC, -1J*delta, 1)
       
          OL = np.abs(tn.overlap(l0,G0, l,G))[0][0]
          rt = -np.log(np.abs(OL)**2)/L
          f1.write(str((i+1)*delta)+"\t"+str(OL.real)+"\t"+str(rt.real)+"\n")
          f2.write(str(disc)+"\n")
          f1.flush()
          f2.flush()
          #print((i+1)*delta, disc)
          disc_tot +=disc
    print(' total discar : ',disc_tot)
    f1.close()
    f2.close()






################################
print("Fig.1==> Phys. Rev. B 87, 195104 (2013)")
print("     ")
start = time.time()
l, G = ground_state()
print('\n elapsed time : ', time.time()-start)
print("\n ground_state is found")
print("     ")
print("     ")
start = time.time()
quench(l, G)
print('\n elapsed time : ', time.time()-start)
print(" \n JOB IS DONE -:) ")

