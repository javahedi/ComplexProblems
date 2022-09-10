#from numba import jit
import numpy as np
from scipy import linalg
import time
import matplotlib.pyplot as plt

import sys,os
#lanczos_path = os.path.join(os.getcwd(),"lanczos")
#sys.path.insert(0,lanczos_path)
#from lanczos import lanczos_full,lanczos_iter,lin_comb_Q_T,expm_lanczos,FTLM_static_iteration

#os.environ['KMP_DUPLICATE_LIB_OK']='True' # uncomment this line if omp error occurs on OSX for python 3
#os.environ['OMP_NUM_THREADS']='4' # set number of OpenMP threads to run in parallel
os.environ['MKL_NUM_THREADS']='4' # set number of MKL threads to run in parallel

import gen_operator as GO
import gen_op_total as GOT
import gen_diagprojector as GD
import gen_nn_int as Gint
#import gen_diagonal_ME as GME
import fast_loop
        



################################################
################################################
################################################
################################################
################################################
def RK4(fk,ham,dt):
    out = fk
    for k in range(4):
        fk *= -1j*dt*ham/(k+1)
        out +=fk
    return out

################################################
################################################
################################################
################################################
################################################

L = 8 #int(np.loadtxt('ChainLength.txt'))
p = 3 #int(np.loadtxt('Flavour.txt'))

projection = True
# %% Generate operators
f0, fu, nop = GO.gen_operator(L,p)
Ntot        = GOT.gen_op_total(nop)
print(" Operators are created ", flush=True)
print( ' Hilber dim : ', f0[0].shape, flush=True)



val = L//2    #np.loadtxt('ParticleNumber.txt')
P   = GD.gen_diagprojector(Ntot.diagonal(), val)
print(" Projection operator is created ", flush=True)
print( ' Projected dim : ', P.shape, flush=True)


class transport():

    def __init__(self):
        pass
    def __str__(self):
        pass
    def autocorr_ED():
        pass
    def autocorr_QDT():
        pass
    def mean_square_distance():
        pass
    def weighted_correlation():
        pass



################################################
################################################
################################################
################################################
################################################

W             = 1.0                                              # np.loadtxt('Wlist')   
g             = 0.0                                              # disorder strenght
g_list        = [g]*L                                            # np.loadtxt('Glist')
m_list        = np.random.uniform(-W, W, L)                      # generate  random on-site
Hint          = Gint.gen_nn_int(fu,nop,g_list, m_list,bc='obc')  # creat hamiltonian
Hint          = P @ Hint @ P.T                                   # projecting Hmailtonian
energy, evecs = linalg.eigh(Hint.todense())                      # full diagonaization of H


        

"""
####################################################
####################################################
####################################################
tmax = 100
timeList = np.linspace(0.,tmax,2000,endpoint=False)
dt = timeList[1] - timeList[0]


Nl = nop[L//2]
dim = f0[0].shape[0]
mu,sigma=0,1
a = np.random.normal(mu,sigma,size=dim)
b = np.random.normal(mu,sigma,size=dim)
phi = (a  + 1j*b)
#phi = np.reshape(phi,(phi.shape[0],1))

if projection:
   phi = P @ phi
   Nl = P @ Nl @ P.T

phi /=np.linalg.norm(phi)
phi = np.sqrt(Nl) @ phi/np.sqrt(np.linalg.norm(phi))
phi2= np.sqrt(Nl) @ phi/np.sqrt(np.linalg.norm(phi))

autocorr_QDT  = []
autocorr_QDT2 = []

start = time.time()
for id, t in enumerate(timeList):
    ##  rk4
    phi = RK4(phi,Hint,dt)
    autocorr_QDT.append((phi.conj().T @ Nl @ phi))
    

    '''
    ##  Krylov
    start = time.time()
    # left
    E, V, Q_T = lanczos_iter(hamiltonian,phi2,20)
    phi2 = expm_lanczos(E,V,Q_T,a=-1j*dt)
    autocorr_QDT2.append( (phi2.conj().T @ Nl @ phi2)) 
    #print(id, 'time is elapse :', time.time()-start)
    '''

print(' autocorre with DQT takes :', time.time()-start,'  seconds')

#####################################################
#####################################################
#####################################################
#####################################################
"""
start = time.time()
Nl = nop[L//2]
Nl = P @ Nl @ P.T

tmax = 100
timeList = np.linspace(0.,tmax,201,endpoint=True)

cc = evecs.conj().T @ Nl @ evecs
#Cython fast loop 
autocorr_ED = fast_loop.TimeDomian(timeList,energy,cc)
autocorr_ED /=energy.shape[0]
print(' autoco with ED takes [CYTHON] :', time.time()-start,'  seconds')

np.savetxt('Autocorr_sample.txt',autocorr_ED.real)
"""
#####################################################
#####################################################
#####################################################
#####################################################
f1 = open('autocorrelation.txt','w')
for i in range(len(timeList)):
    f1.write('%0.8f     %0.8f      %0.8f\n'%(timeList[i],autocorr_QDT[i].real,autocorr_ED[i].real))
f1.close()
"""
####################################################
####################################################
####################################################
#Christoph Karrasch
#arXiv:1702.04349v2
'''
x2 = 0
for i in range(L):
    Ni =  P @ ( fu[i].T @ fu[i])@ P.T # -0.5*f0[n] ) @ P.T
    ci =  V @ Ni @ V.T #np.einsum('ij,jk,kl->il',V,Ni,V.T)
    for j in range(L):
        Nj =  P @ ( fu[j].T @ fu[j])@ P.T # -0.5*f0[m] ) @ P.T
        cj = V @ Nj @ V.T #np.einsum('ij,jk,kl->il',V,Nj,V.T)
        
        Cij = 0
        for n,en in enumerate(E):
            for m,em in enumerate(E):
                #cij = np.sum(V[n,:].T @ Ni @ V[m,:]) * np.sum(V[m,:].T @ Nj @ V[n,:])
                Cij +=  np.exp(1j*(en-em)*t) * ci[n,m] * cj[m,n]
                #print(cn[n,m] * cm[m,n]-cij)
        x2 += (i-j)**2 * Cij
                
x2 /=L
x2 /=E.shape[0]

print(r,'  time :', time.time() - start)

np.savetxt('x2_'+str(r)+'.txt',x2.real)
MSD_samples +=x2.real
'''
start = time.time()

#Time = np.linspace(0.01,10,101,endpoint=False)
Time = np.logspace(0,1.7,101,base=10)


N0 =  P @ nop[0] @ P.T
c0 =  evecs.conj().T @ N0 @ evecs
x2 = 0
   
for r in range(1,L):
    Nr =  P @ nop[r] @ P.T
    cr =  evecs.conj().T @ Nr @ evecs
    Cr0_0 = np.sum(evecs.conj().T @ Nr @ N0 @ evecs)
    
    #Cr0_t = 0
    #for i,ei in enumerate(E):
    #    for j,ej in enumerate(E):
    #        Cr0_t +=  np.exp(1j*(ej-ei)*Time) * cr[i,j] * c0[j,i]

    Cr0_t = fast_loop.TimeDomian2(Time,energy,cr,c0)
    x2 += r**2 * (Cr0_t-Cr0_0*np.ones_like(Time))
            

x2 /=energy.shape[0]
np.savetxt('MSD_sample.txt',x2.real)

print(' MSD with ED takes [CYTHON] :', time.time()-start,'  seconds')

print(" JOB IS FINISHED :-) ")
