import numpy as np
from copy import deepcopy
import tn
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
    
    
    ######## Define the simulation parameter ######################
    #OBC = 0; IBC =1
    BC = 0
    d  = 4
    chiMax = 10
    delta = 0.01
    ######## Define the model parameter ###########################
    J = 1.
    g = 0.0
    L = 100
    Jxy1 = (L-1)*[J]
    Jz1  = (L-1)*[J]
    Jxy2 = (L-2)*[J]
    Jz2  = (L-2)*[J]
    hz = L*[0.]
    hx = L*[g]
    N = int(5./delta)
    ########   Initial MPS     ####################################
    l, G, chi= tn.initialMat(L//2,d,chiMax,BC)
    l, G, norm = tn.normalize_left (l,  G)
    
    expH_bond,H_bond = tn.mixing_hamiltinian_expH_bond(Jxy1,Jz1,Jxy2,Jz2,hz,hx,L,delta)
    for i in tqdm.tqdm(range(N)):
        #print('%i/%i'%(i,N))
        tn.Suzuki_Trotter_1ftOrder_two_sites(l,G,expH_bond,False)
    sum = 0
    for pos in range(L//2-1):
        sum += tn.expc_two_site(l,G,H_bond[pos],pos).real
        
    print("\n Ground state energy per spin  : ", sum/L)
    
        
        



ground_state()

