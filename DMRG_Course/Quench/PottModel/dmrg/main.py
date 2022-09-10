import numpy as np
from copy import deepcopy
import tn
import tqdm
import time
import os, sys
os.environ['MKL_NUM_THREADS']='2' # set number of MKL threads to run in parallel
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
    BC = 0
    d  = 3
    chimax = 40
    delta = 0.01
    ######## Define the model parameter ###########################
    L = 20
    J = L*[-0.0]
    f = L*[-1.0]
    
    N = int(5./delta)
    ########   Initial MPS     ####################################
    l, G, chi= tn.initialMat(L,d,chimax,BC)
    l, G, norm = tn.normalize_left (l,  G)
    
    expH_bond,H_bond = tn.pott_expH_bond(J,f,L,delta)
    for i in tqdm.tqdm(range(N)):
        #print('%i/%i'%(i,N))
        tn.Suzuki_Trotter_1ftOrder_two_sites(l,G,expH_bond,infinite)
        
    sum = 0.
    for pos in range(L-1):
        sum += tn.expc_two_site(l,G,H_bond,pos).real
    
    print(sum/L,"(DMRG with dt : ", delta,")")
    return l, G
        
        



def quench(l, G):
     l0, G0 = deepcopy(l), deepcopy(G)
     infinite = False
     ######## Define the model parameter ###########################
     L = 20
     chimax = 60
     J = L*[-1.0]
     f = L*[-0.0]
     T = 4
     delta = 0.01
     N = int(T/delta)
     ######### Return function  ##########
     expH_bond,H_bond = tn.pott_expH_bond(J,f,L,-1j*delta)
     rt = []
     t  = []
     f1 = open('Pott_L_'+str(L)+'_chi_'+str(chimax)+'.txt','w')
     for i in tqdm.tqdm(range(N)):
         tn.Suzuki_Trotter_1ftOrder_two_sites(l,G,expH_bond,infinite)
         t.append((i+1)*delta)
         rt.append(-np.log(np.abs(tn.overlap_v2(l0,G0, l,G))**2)/L)
         f1.write('%f    %f     %f\n'%(t[-1],np.abs(tn.overlap_v2(l0,G0, l,G)), rt[-1].real))
         
     f1.close()
     
     

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


