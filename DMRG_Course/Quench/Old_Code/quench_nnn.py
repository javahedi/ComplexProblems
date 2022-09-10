import numpy as np
from copy import deepcopy
import tn
from ed_ising import ising_gs
import tqdm
import time
import os, sys
#os.environ['MKL_NUM_THREADS']='2' # set number of MKL threads to run in parallel
#os.environ['OMP_NUM_THREADS']='1' # set number of OpenMP threads to run in parallel


import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


#Set global matplotlib parameters in script or in /home/$USER/.matplotlib/matplotlibrc
plt.rcParams['axes.linewidth'] = 1.
plt.rcParams['xtick.major.size'] = 4
plt.rcParams['xtick.minor.size'] = 2
plt.rcParams['ytick.major.size'] = 4
plt.rcParams['ytick.minor.size'] = 2
plt.rcParams['xtick.bottom']       = True
plt.rcParams['xtick.top']          = True
plt.rcParams['ytick.left']         = True
plt.rcParams['ytick.right']        = True
plt.rcParams['xtick.direction']    = 'in'
plt.rcParams['ytick.direction']    = 'in'
plt.rcParams['font.family']        = 'sans-serif'

###############################################


def ground_state():
    
    infinite=False
    ######## Define the simulation parameter ######################
    #OBC = 0; IBC =1
    BC = 0
    d  = 2
    chiMax = 100
    delta = 0.01
    ######## Define the model parameter ###########################
    J = -1.
    g = 1.3
    if infinite:
       L = 2
       Jxy1 = L*[0.]
       Jz1  = L*[J]
       Jxy2 = L*[0.]
       Jz2  = L*[0.]
    else:
       L = 15
       Jxy1 = (L-1)*[0.]
       Jz1  = (L-1)*[J]
       Jxy2 = (L-1)*[0.]
       Jz2  = (L-1)*[0.]
    hz = L*[0.]
    hx = L*[g]
    N = int(5./delta)
    ########   Initial MPS     ####################################
    l, G, chi= tn.initialMat(L,d,chiMax,BC)
    l, G, norm = tn.normalize_left (l,  G)
    
    expH_bond,H_bond = tn.nnn_hamiltinian_expH_bond(Jxy1,Jz1,Jxy2,Jz2,hz,hx,L,delta)
    for i in tqdm.tqdm(range(N)):
        #print('%i/%i'%(i,N))
        tn.Suzuki_Trotter_1ftOrder_three_sites(l,G,expH_bond,infinite)
    return l, G
        
        



def quench(l, G):
     l0, G0 = deepcopy(l), deepcopy(G)
     
     infinite=False
     
     ######## Define the model parameter ###########################
     J = -1.
     g = 0.2
     if infinite:
       L = 2
       Jxy1 = L*[0.]
       Jz1  = L*[J]
       Jxy2 = L*[0.]
       Jz2  = L*[0.6]
     else:
       L = 15
       Jxy1 = (L-1)*[0.]
       Jz1  = (L-1)*[J]
       Jxy2 = (L-2)*[0.]
       Jz2  = (L-2)*[0.6]
     hz = L*[0.]
     hx = L*[g]
     T = 5
     delta = 0.01
     N = int(T/delta)
     ######### Return function  ##########
     expH_bond,H_bond = tn.nnn_hamiltinian_expH_bond(Jxy1,Jz1,Jxy2,Jz2,hz,hx,L,delta)
     rt = []
     t  = []
     f1 = open('quench_PRB_L_'+str(L)+'_nnn.txt','w')
     for i in tqdm.tqdm(range(N)):
         tn.Suzuki_Trotter_1ftOrder_three_sites(l,G,expH_bond,infinite)
         t.append((i+1)*delta)
         rt.append(-np.log(np.abs(tn.overlap_v2(l0,G0, l,G))**2)/L)
         f1.write('%f    %f     %f\n'%(t[-1],np.abs(tn.overlap_v2(l0,G0, l,G)), rt[-1].real))
         
     f1.close()
     
     

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


