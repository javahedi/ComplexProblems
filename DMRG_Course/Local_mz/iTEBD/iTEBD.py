import numpy as np 
from scipy.linalg import expm
from svd_robust import svd
import tqdm
from scipy import integrate


##############################
#        Creating MPS
###############################
def initialMat(L,d,chiMax,bc='OBC'):
    G = []  # Gamma
    l = []  # Lambda
    chi = np.zeros(L+1, dtype=np.int32)
    
    if bc=='OBC':
        for pos in range(0,L//2+1):
            chi[pos] = np.minimum( chiMax , d**pos     )
        for pos in range(L//2+1,L+1):
            chi[pos] = np.minimum( chiMax , d**(L-pos) )
    else:
        chi[:] = chiMax
 
    for pos in range(L):
        chi0 = int( chi[ pos    ] )
        chi1 = int( chi[ pos+1  ] )
        G.append( np.zeros( (d,chi0,chi1), np.complex128 ) )
        l.append( np.zeros( (  chi0     )                ) )
    
    chi0 = int( chi[ L ] )
    l.append(np.zeros( (  chi0     )       ))
    ########################################
    #       Initial state
    ########################################
    for pos in range(L):
        G[pos][pos%2,0,0] = 1.0 # example: Neel state
        #G[pos][0   ,0,0] = 1.0 # example: FM state
        l[pos][   :     ] = 1.0
    l[L][   :     ] = 1.0
    return G, l
    
    
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
        #if i%2==1:
        #    G[-1][0,0,0]=1
        #else:
        #    G[-1][1,0,0]=1
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
    sx = np.array([[0., 1.],[1. , 0.]])/2.
    sy = np.array([[0.,-1j],[1j , 0.]])/2.
    sz = np.array([[1., 0.],[0. ,-1.]])/2.
    d = 2

    U_bond = []
    H_bond = []
    for i in range(L):
        H = J * np.real(np.kron(sx,sx) + np.kron(sy,sy)) + D*np.kron(sz,sz)
        H += h*  np.kron(sz,np.eye(2))
        H_bond.append(np.reshape(H,(d,d,d,d)))
        U_bond.append(np.reshape(expm(-delta*H),(d,d,d,d)))

    return U_bond,H_bond

########################################
#       Creating Matrices A, B
########################################
def lG_AB(l, G):
    L    = len(G)
    A = []
    B = []
    for pos in range(L):
        temp = np.tensordot( np.diag(l[pos]), G[pos], axes=(1,1) )
        temp = np.transpose( temp, (1,0,2) )
        A.append( temp )
        
        temp = np.tensordot( G[pos], np.diag(l[pos+1]), axes=(2,0) )
        B.append( temp )
    return A, B
    
########################################
#    Left-Normalization for Vidal MPS
########################################
def normalize_left(l, G):
    A, B = lG_AB(l, G)
    AB   = A
    L    = len(G)
    d    = B[0].shape[0]
    
    for pos in range(L):
        chi0 = AB[pos].shape[1]
        chi1 = AB[pos].shape[2]
        if pos == 0 :
           M        = np.reshape   ( B[pos], (d*chi0,chi1)   )
           U, S, V  = svd          ( M     , full_matrices=0 )
           U        = np.reshape   ( U     , (d,chi0,chi1)   )
           G[pos]   = U
        else        :
           SV       = np.tensordot ( np.diag(S)                 , V              , axes=(1,0) )
           M        = np.tensordot ( SV                         , B[pos]         , axes=(1,1) )
           M        = np.transpose ( M                          , (1,0,2)                     )
           M        = np.reshape   ( M                          , (d*chi0,chi1)               )
           U, S, V  = svd          ( M                          , full_matrices=0             )
           U        = np.reshape   ( U                          , (d,chi0,chi1)               )
           linvU    = np.tensordot ( np.diag(inverse(l[pos]))   , U              , axes=(1,1) )
           linvU    = np.transpose ( linvU                      , (1 ,0 ,2)                   )
           G[pos]   = linvU
           
        
        if pos < L-1 : l[pos+1] = S
        else         : norm = (np.tensordot( np.diag(S), V, axes=(1,0))[0,0]).real
    return l, G, norm

########################################
#    Right-Normalization for Vidal MPS
########################################
def normalize_right(l, G):
    A, B = lG_AB(l, G)
    AB   = B
    L    = len(G)
    d    = A[0].shape[0]
    
    for pos in range(L-1,-1,-1):
        chi0 = AB[pos].shape[1]#chi[pos]
        chi1 = AB[pos].shape[2]#chi[pos+1]
        if pos == L-1 :
           M        = np.transpose ( A[pos], (1,0,2)         )
           M        = np.reshape   ( M     , (chi0,d*chi1)   )
           U, S, V  = svd          ( M     , full_matrices=0 )
           V        = np.reshape   ( V     , (chi0,d,chi1)   )
           V        = np.transpose ( V     , (1,0,2)         )
           G[pos]   = V
        else        :
           US       = np.tensordot ( U     , np.diag(S)                   , axes=(1,0) )
           M        = np.tensordot ( A[pos], US                           , axes=(2,0) )
           M        = np.transpose ( M     , (1,0,2)                                   )
           M        = np.reshape   ( M     , (chi0,d*chi1)                             )
           U, S, V  = svd          ( M     , full_matrices=0                           )
           V        = np.reshape   ( V     , (chi0,d,chi1)                             )
           V        = np.transpose ( V     , (1,0,2)                                   )
           Vlinv    = np.tensordot ( V     , np.diag(inverse(l[pos+1]))   , axes=(2,0) )
           G[pos]   = Vlinv

        if pos > 0 :  l[pos] = S
        else       :  norm = (np.tensordot( U, np.diag(S), axes=(1,0))[0,0]).real
    return l, G, norm
               


########################################
#           MPS overlap
########################################
def overlap(l1, G1, l2, G2):
    L = len(G1)
    A1, B1 = lG_AB(l1, G1)
    A2, B2 = lG_AB(l2, G2)
        
    pos = 0
    OL_now =  np.tensordot (np.conj(A1[pos]),  A2[pos]          , axes=([0,1],[0,1]) )
    for pos in range(1,L):
        temp = np.tensordot(OL_now          ,  np.conj(A1[pos]) , axes=(0,1)         )
        temp = np.tensordot(temp            ,  A2[pos]          , axes=([1,0],[0,1]) )
        OL_now = temp
        
    return OL_now
        

###########################################################
#     iTEBD   sweep
##########################################################
def itebd(G,l,expH,chimax):
    d = G[0].shape[0]
    L = len(G)
    for pos in [0,1]:
        ia   = np.mod(pos  ,L)
        ib   = np.mod(pos+1,L) ## ib==pos
        chi1 = G[ia].shape[1]
        #chi2 = G[ia].shape[2]
        chi3 = G[ib].shape[2]
 
        # construct psi
        psi = np.tensordot(np.diag(l[ib]), G[ia]          , axes=(1,1)) 
        psi = np.tensordot(psi           , np.diag(l[ia]) , axes=(2,0))
        psi = np.tensordot(psi           , G[ib]          , axes=(2,1)) 
        psi = np.tensordot(psi           , np.diag(l[ib]) , axes=(3,0))

        # apply expH
        psi = np.tensordot(psi , expH[ib], axes=([1,2],[0,1]) )
        
        # SVD 
        psi     = np.transpose(psi, (2,0,3,1)       )
        psi     = np.reshape  (psi, (d*chi1,d*chi3) )
        U, S, V = svd(psi)
        V = V.T
        #V = V.T  V_(chi2)(d*chi3)==>V_(d*chi3)(chi2)
        chi2    = np.max([np.sum(S>10**(-10)),chimax])

        # truncate
        l[ia] = S[0:chi2]/np.sqrt(np.sum(S[0:chi2]**2) )
      
        U = np.reshape  (U[:,0:chi2]         , (d,chi1,chi2)             )
        U = np.tensordot(np.diag(l[ib]**(-1)), U             ,axes=(1,1) )
        U = np.transpose(U                   , (1,0,2)                   )
        G[ia] = U

        V = np.reshape  (V[:,0:chi2] , (d,chi3,chi2)                        )
        V = np.transpose(V           , (0,2,1)                              )
        V = np.tensordot(V           ,  np.diag(l[ib]**(-1)) ,   axes=(2,0) )
        G[ib] = V



########################################
#    two-sites expection (nearest-neighbor)
########################################
def two_expectation2(G, l, OP):
    " Expectation value for a bound operators "
    E  = []
    L = len(G)
    for pos in range(L):
        ia  = np.mod(pos  ,L)
        ib  = np.mod(pos+1,L)
        Psi = np.tensordot( np.diag(l[ib])  ,  G[ia]          , axes=(1,1)                 )
        Psi = np.tensordot( Psi             ,  np.diag(l[ia]) , axes=(2,0)                 )
        Psi = np.tensordot( Psi             ,  G[ib]          , axes=(2,1)                 )
        Psi = np.tensordot( Psi             ,  np.diag(l[ib]) , axes=(3,0)                 )
        Psi = np.tensordot( np.conj(Psi)    ,  Psi            , axes=([0,3],[0,3])         )
        E.append(np.squeeze(np.tensordot(Psi,  OP[ib]         , axes=([0,1,2,3],[0,1,2,3]) )).item())
        
    return E
    
def two_expectation(G, l ,OP):
    " Expectation value for a bound operators "
    E  = []
    L = len(G)
    for pos in range(L):
        ia     = np.mod(pos  ,L)
        ib     = np.mod(pos+1,L)
        psi    = np.tensordot(np.diag(l[ib]), G[ia]          , axes=(1,1)         )
        psi    = np.tensordot(psi           , np.diag(l[ia]) , axes=(2,0)         )
        psi    = np.tensordot(psi           , G[ib]          , axes=(2,1)         )
        psi    = np.tensordot(psi           , np.diag(l[ib]) , axes=(3,0)         )
        psi_OP = np.tensordot(psi           , OP[ib]     , axes=([1,2],[0,1]) )
        psi_OP = psi_OP.conj()
        E.append(np.squeeze(np.tensordot(psi_OP,psi,axes=([0,1,2,3],[0,3,1,2]))).item())
    return E


########################################
#           one-site expection
########################################
def one_expectation(G, l, OP):
    " Expectation value for a site operator "
    E  = []
    L = len(G)
    for pos in range(L):
        ia  = np.mod(pos  ,L)
        ib  = np.mod(pos+1,L) # ib==pos
        print(np.diag(l[ib]).shape, G[ib] .shape)
        Psi = np.tensordot( np.diag(l[ib])   ,  G[ib]          , axes=(1,1)         )
        Psi = np.tensordot( Psi              ,  np.diag(l[ia]) , axes=(2,0)         )
        Psi = np.tensordot( np.conj(Psi)     ,  Psi            , axes=([0,2],[0,2]) )
        E.append(np.squeeze(np.tensordot(Psi ,  OP[ib]         , axes=([0,1],[0,1]) )).item())
    return E
    
###############################################
def local_dynamic():
    ######## Define the simulation parameter ######################
    L      = 2
    d      = 2
    chimax = 10
    sz     = np.array([[1.,0.],[0.,-1.]])
    ######## Define the model parameter ###########################
    J     = 1.0
    D     = 0.5
    g     = 0.0
    delta  = 0.1
    tfin   = 10.0
    N      = int(tfin/delta)
    ######### initial MPS state ##########
    #G,l = init_fm_mps(L)
    G, l = initialMat(L,d,chimax,bc='IBC')
    expH,H = heisenberg(J,D,g,L,1j*delta)
    
    f1 = open('mz_chimax'+str(chimax)+'_IBC.txt','w')
    #for i in tqdm.tqdm(range(N)):
    for i in range(N):
        itebd(G,l,expH,chimax)
        m = one_expectation(G,l,[sz,sz])[0]

        
        f1.write(str((i+1)*delta)+"\t"+str(m.real)+"\n")
    f1.flush()
          
    f1.close()
    

local_dynamic()

