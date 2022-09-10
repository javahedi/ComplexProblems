import numpy as np 
from scipy.linalg import expm
from svd_robust import svd
import tqdm
from scipy import integrate



def init_fm_mps(L):
    """ Returns FM  MPS"""
    d = 2
    G = []
    l = []
    for i in range(L):
        G.append(np.zeros([2,1,1]))
        G[-1][0,0,0]=1
        l.append(np.ones([1]))
    l.append(np.ones([1]))
    return G,l
    
def init_af_mps(L):
    """ Returns AF  MPS"""
    d = 2
    G = []
    l = []
    for i in range(L):
        G.append(np.zeros([2,1,1]))
        G[-1][i%2,0,0]=1
        l.append(np.ones([1]))
    l.append(np.ones([1]))
    return G,l
    
def ising(g,J,L,delta):
    """ Returns bond hamiltonian and bond time-evolution"""
    sx = np.array([[0.,1.],[1.,0.]])
    sy = np.array([[0.,-1j],[1j,0.]])
    sz = np.array([[1.,0.],[0.,-1.]])
    d = 2

    U_bond = []
    H_bond = []
    for i in range(L):
        H = J*np.kron(sz,sz) + g*np.kron(sx,np.eye(2))
        H_bond.append(np.reshape(H,(d,d,d,d)))
        U_bond.append(np.reshape(expm(-delta*H),(d,d,d,d)))

    return U_bond,H_bond

def heisenberg(J,D,h,L,delta):
    """ Returns bond hamiltonian and bond time-evolution"""
    sx = np.array([[0.,1.],[1.,0.]])/2.
    sy = np.array([[0.,-1j],[1j,0.]])/2.
    sz = np.array([[1.,0.],[0.,-1.]])/2.
    d = 2

    U_bond = []
    H_bond = []
    for i in range(L):
        H = J * np.real(np.kron(sx,sx) + np.kron(sy,sy)) + D*np.kron(sz,sz)
        H += h*  np.kron(sz,np.eye(2))
        H_bond.append(np.reshape(H,(d,d,d,d)))
        U_bond.append(np.reshape(expm(-delta*H),(d,d,d,d)))

    return U_bond,H_bond


def itebd(G,l,expH,chimax):
    d = G[0].shape[0]
    for ibound in [0,1]:
        ia   = np.mod(ibound,2)
        ib   = np.mod(ibound+1,2)
        chi1 = G[ia].shape[1]
        chi3 = G[ib].shape[2]
 
        # construct psi
        psi = np.tensordot(np.diag(l[ib]), G[ia]          , axes=(1,1)) 
        psi = np.tensordot(psi           , np.diag(l[ia]) , axes=(2,0))
        psi = np.tensordot(psi           , G[ib]          , axes=(2,1)) 
        psi = np.tensordot(psi           , np.diag(l[ib]) , axes=(3,0))

        # apply expH
        psi = np.tensordot(psi , expH[ibound], axes=([1,2],[0,1] ))
        
        # SVD 
        psi     = np.transpose(psi, (2,0,3,1)       )
        psi     = np.reshape  (psi, (d*chi1,d*chi3) )
        X, Y, Z = svd(psi)
        Z = Z.T
        #Z = Z.T  Z_(chi2)(d*chi3)==>Z_(d*chi3)(chi2)
        chi2    = np.min([np.sum(Y>10**(-10)),chimax])

        # truncate
        l[ia] = Y[0:chi2]/np.sqrt(np.sum(Y[0:chi2]**2) )
      
        X = np.reshape  (X[:,0:chi2]         , (d,chi1,chi2)             )
        X = np.tensordot(np.diag(l[ib]**(-1)), X             ,axes=(1,1) )
        X = np.transpose(X                   , (1,0,2)                   )
        G[ia] = X

        Z = np.reshape  (Z[:,0:chi2] , (d,chi3,chi2)                        )
        Z = np.transpose(Z           , (0,2,1)                              )
        Z = np.tensordot(Z           ,  np.diag(l[ib]**(-1)) ,   axes=(2,0) )
        G[ib] = Z


########################################
#           one-site expection
########################################
def expvalue1(G ,l ,OP):
    " Expectation value for a site operator "
    E  = []
    L = len(G)
    for pos in range(L):
        ia     = np.mod(pos,L)
        ib     = np.mod(pos+1,L)  # ib==pos
        Psi = np.tensordot( np.diag(l[ib])   ,  G[ib]          , axes=(1,1)         )
        Psi = np.tensordot( Psi              ,  np.diag(l[ia]) , axes=(2,0)         )
        Psi = np.tensordot( np.conj(Psi)     ,  Psi            , axes=([0,2],[0,2]) )
        E.append(np.squeeze(np.tensordot(Psi ,  OP[ib]         , axes=([0,1],[0,1]) )).item())
        
    return E
########################################
#    two-sites expection (nearest-neighbor)
########################################
def expvalue2(G,l,OP):
    " Expectation value for a bound operators "
    E  = []
    L = len(G)
    for pos in range(L):
        ia     = np.mod(pos,L)
        ib     = np.mod(pos+1,L)
        Psi = np.tensordot( np.diag(l[ib])  ,  G[ia]          , axes=(1,1)                 )
        Psi = np.tensordot( Psi             ,  np.diag(l[ia]) , axes=(2,0)                 )
        Psi = np.tensordot( Psi             ,  G[ib]          , axes=(2,1)                 )
        Psi = np.tensordot( Psi             ,  np.diag(l[ib]) , axes=(3,0)                 )
        Psi = np.tensordot( np.conj(Psi)    ,  Psi            , axes=([0,3],[0,3])         )
        E.append(np.squeeze(np.tensordot(Psi,  OP[ib]         , axes=([0,1,2,3],[0,1,2,3]) )).item())
        
    return E
    
def bond_expectation(G ,l ,OP):
    " Expectation value for a bound operators "
    E  = []
    L = len(G)
    for ibound in range(L):
        ia     = np.mod(ibound,L)
        ib     = np.mod(ibound+1,L)
        psi    = np.tensordot(np.diag(l[ib]), G[ia]          , axes=(1,1)         )
        psi    = np.tensordot(psi           , np.diag(l[ia]) , axes=(2,0)         )
        psi    = np.tensordot(psi           , G[ib]          , axes=(2,1)         )
        psi    = np.tensordot(psi           , np.diag(l[ib]) , axes=(3,0)         )
        psi_OP = np.tensordot(psi           , OP[ibound]     , axes=([1,2],[0,1]) )
        psi_OP = psi_OP.conj()
        E.append(np.squeeze(np.tensordot(psi_OP,psi,axes=([0,1,2,3],[0,3,1,2]))).item())
    return E



######## Define the model and simulation parameters ######################
chimax = 20
d      = 2
L      = 2
g      = 0.5
J      = -1.0
D      = 0.0
sz = np.array([[1.,0.],[0.,-1.]])


for delta in [0.1,0.01,0.001]:
    N = int(20./delta)
    G,l = init_af_mps(L)
    U_bond,H_bond = heisenberg(J,D,g,L,delta)
    #U_bond,H_bond = ising(g,J,L,delta)
    #for i in tqdm.tqdm(range(N)):
    for i in range(N):
        itebd(G,l,U_bond,chimax)
    E = np.mean(expvalue2(G,l,H_bond))
    #E2 = np.mean(bond_expectation(G,l,H_bond))
        
    #m = np.mean(expvalue1(G,l,[sz,sz]))
    #S = np.mean(entanglement_entropy(s))
       

    print(E,"(iTEBD, delta =" ,delta, ")")


# Ising
#f = lambda k,g : -2*np.sqrt(1+g**2-2*g*np.cos(k))/np.pi/2.
#E0_exact = integrate.quad(f, 0, np.pi, args=(g,))[0]
# XX
E0_exact = -1.0/np.pi
print(E0_exact,"(EXACT)")
