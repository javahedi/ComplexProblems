import numpy as np
from copy import deepcopy
import tn
from ed_ising import ising_gs
from scipy import integrate
import tqdm


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

#OBC = 0; IBC =1
#L = 20
#d = 2
#chiMax = 20
#BC = OBC

#l, G, chi= tn.initialMat(L,d,chiMax,0)
#lini = deepcopy(l)
#Gini = deepcopy(G)

#l, G, norm = tn.normalize_right (l,  G)


#normini = np.sqrt(np.abs(tn.overlap(lini, Gini, lini, Gini)))
#print(normini)
#print((tn.overlap(l, G, l, G)/normini).real)
#print('%$'*20)
#for pos in range(L+1):
    #Psi = np.tensordot( np.diag(l[pos]) , G[pos] , axes=(1,1) )
    #print(l[pos])






def clean_quench():
     ######## Define the simulation parameter ######################
     OBC = 0; IBC =1
     L = 20
     d = 2
     chiMax = 300
     BC = OBC
     f1 = open('timeev_L20_dt01_chi300.txt','w')
     ######## Define the model parameter ###########################
     Jxy = (L-1)*[1.]
     Jz  = (L-1)*[0.5]
     hz  = L*[0.]
     hx  = L*[0.]
     delta = 0.1
     T   = 10
     N   = int(T/delta)
     sz  = np.array([[1.,0.],[0.,-1.]])/2.
     
     ######### initial MPS state ##########
     l, G, chi= tn.initialMat(L, d, chiMax, BC)
     #l, G, norm = tn.normalize_left (l,  G)
     #l, G, norm = tn.normalize_right(l,  G)

       
     
     G0 = G
     l0 = l
     
     expH_bond,H_bond = tn.nn_hamiltinian_expH_bond(Jxy,Jz,hz,hx,L,delta*1j)
     S = []
     m = []
     LE = []
     bond = []
     t = []
     for i in tqdm.tqdm(range(N)):
         tn.Suzuki_Trotter_1ftOrder_two_sites(l,G,expH_bond,infinite=False)
         t.append((i+1)*delta)
         S.append(tn.entanglement_entropy(l)[L//2])
         m.append(tn.expc_one_site(l,G,sz,L//2).real)
         LE.append(np.abs(np.sum(tn.loschmidt_echo_v2(l0,G0, l,G)))**2)
         
         f1.write('%f    %f\n'%(t[-1],2*m[-1].real))
         
         #print(i,LE[-1])
         
        
         

     f1.close()
     ######## Plot ##########################
     
     plt.figure()
     plt.plot(t,m)
     plt.ylim([np.min(m)-0.1,np.max(m)+0.1])
     plt.ylabel(r'$m^z$')
     plt.xlabel(r'$t$')
     plt.show()
     
     
     plt.figure()
     plt.plot(t,S)
     plt.ylim([np.min(S)-0.1,np.max(S)+0.1])
     plt.ylabel(r'$S(t)$')
     plt.xlabel(r'$t$')
     plt.show()
     
     

        
        




def ground_state():

    infinite=False
    
     
    ######## Define the simulation parameter ######################
    OBC = 0; IBC =1
    d = 2
    chiMax = 20
    BC = OBC
    ######## Define the model parameter ###########################
    J = 1.
    g = 1.
    if infinite:
        L = 2
        Jxy = L*[0.]
        Jz =  L*[J]
    else:
        L = 8
        Jxy = (L-1)*[0.]
        Jz  = (L-1)*[J]
    
    
    hz = L*[0.]
    hx = L*[g]
    
    sz = np.array([[1.,0.],[0.,-1.]])
    szsz = tn.op_tensor(sz,sz)
    

    l, G, chi  = tn.initialMat (L, d, chiMax, BC)
    l, G, norm = tn.normalize_left (l,  G)

    for delta in [0.1,0.01,0.001]:
        N = int(5./delta)
        
        expH_bond,H_bond = tn.nn_hamiltinian_expH_bond(Jxy,Jz,hz,hx,L,delta)

        for i in range(N):
            tn.Suzuki_Trotter_1ftOrder_two_sites(l,G,expH_bond)
        sum = 0
        for pos in range(L-1):
            sum += tn.expc_two_site(l,G,H_bond[pos],pos).real
            #sum += tn.expc_two_site(l,G,szsz,pos).real
        print(sum,"(DMRG with dt : ", delta,")")

    if infinite:
        f = lambda k,g : -2*np.sqrt(1+g**2-2*g*np.cos(k))/np.pi/2.
        E0_exact = integrate.quad(f, 0, np.pi, args=(g,))[0]
        print(E0_exact*2,"(EXACT)")
    else:
        print(ising_gs(g,J,L),"(ED)")
       
    




#clean_quench()
print("    ")
print("Ground state:  using the  imaginary time evolution approach")
ground_state()

