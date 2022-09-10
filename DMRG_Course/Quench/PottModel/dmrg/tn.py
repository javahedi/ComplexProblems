import numpy as np
from scipy import linalg
from svd_robust import svd
from scipy.linalg import expm, sinm, cosm
import sys,os

##############################
#        Creating MPS
###############################
def initialMat(L,d,chiMax,bc='obc'):
    G = []  # Gamma
    l = []  # Lambda
    chi = [0] * (L+1)

    for pos in range(0,L//2+1):
        chi[pos] = np.minimum( chiMax , d**pos     )
    for pos in range(L//2+1,L+1):
        chi[pos] = np.minimum( chiMax , d**(L-pos) )

    for pos in range(L):
        chi0 = int( chi[ pos    ] )
        chi1 = int( chi[ pos+1  ] )
        G.append( np.zeros( (d,chi0,chi1), np.complex128 ) )
        l.append( np.zeros( (  chi0     )                ) )
    
    chi0 = int( chi[ L    ] )
    l.append(np.zeros( (  chi0     )                ))
    ########################################
    #       Initial state
    ########################################
    for pos in range(L):
        #G[pos][pos%2,0,0] = 1.0 # example: Neel state
        G[pos][0   ,0,0] = 1.0 # example: FM state
        l[pos][   :     ] = 1.0
    l[L][   :     ] = 1.0
    return l, G, chi


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



###############################################
#       Local Hamiltonian nearest neighbour
###############################################
def nn_hamiltinian_expH_bond(Jxy,Jz,hz,hx,L,delta):
    """ Returns bond hamiltonian and bond time-evolution"""
    sx = np.array([[0., 1.],[1., 0.]])
    sy = np.array([[0.,-1j],[1j, 0.]])
    sz = np.array([[1., 0.],[0.,-1.]])
    d = 2
    I1 = np.eye(2)

    expH_bond = []
    H_bond = []
   
    for i in range(L-2):
        H =  Jxy[i] * np.real(np.kron(sx,sx) + np.kron(sy,sy))
        H +=  Jz[i] * np.kron(sz,sz)
        H +=  hz[i] * np.kron(sz,I1)
        H +=  hx[i] * np.kron(sx,I1)
        H_bond.append(np.reshape(H,(d,d,d,d)))
        expH_bond.append(np.reshape(expm(-delta*H),(d,d,d,d)))

    i = L-2
    H  = Jxy[i] * np.real(np.kron(sx,sx) + np.kron(sy,sy)) + Jz[i]*np.kron(sz,sz)
    H +=  hz[i] * np.kron(sz,I1) + hz[i+1]*np.kron(I1,sz)
    H +=  hx[i] * np.kron(sx,I1) + hx[i+1]*np.kron(I1,sx)


    H_bond.append(np.reshape(H,(d,d,d,d)))
    expH_bond.append(np.reshape(expm(-delta*H),(d,d,d,d)))
    
    
    return expH_bond,H_bond

###############################################
#       Local Hamiltonian next nearest neighbour
###############################################
def nnn_hamiltinian_expH_bond(Jxy1,Jz1,Jxy2,Jz2,hz,hx,L,delta):
    """ Returns bond hamiltonian and bond time-evolution"""
    sx = np.array([[0., 1.],[1., 0.]])
    sy = np.array([[0.,-1j],[1j, 0.]])
    sz = np.array([[1., 0.],[0.,-1.]])
    d = 2
    I1 = np.eye(2)
    
   

    expH_bond = []
    H_bond    = []
   
    for i in range(L-3):
        H  = Jxy1[i] * np.real(np.kron(np.kron(sx,sx),I1) + np.kron(np.kron(sy,sy),I1))
        H += Jz1[i]  *         np.kron(np.kron(sz,sz),I1)
        
        
        H += Jxy2[i] * np.real(np.kron(np.kron(sx,I1),sx) + np.kron(np.kron(sy,I1),sy))
        H += Jz2[i]  *         np.kron(np.kron(sz,I1),sz)
        
        H += hz[i]   *         np.kron(np.kron(sz,I1),I1)
        H += hx[i]   *         np.kron(np.kron(sx,I1),I1)
        
        H_bond.append(np.reshape(H,(d,d,d,d,d,d)))
        expH_bond.append(np.reshape(expm(-delta*H),(d,d,d,d,d,d)))
    
    i = L-3
    H  = Jxy1[i]   * np.real(np.kron(np.kron(sx,sx),I1) + np.kron(np.kron(sy,sy),I1))
    H += Jz1[i]    *         np.kron(np.kron(sz,sz),I1)
    H += Jxy1[i+1] * np.real(np.kron(np.kron(I1,sx),sx) + np.kron(np.kron(I1,sy),sy))
    H += Jz1[i+1]  *         np.kron(np.kron(I1,sz),sz)
    
    H += Jxy2[i]   * np.real(np.kron(np.kron(sx,I1),sx) + np.kron(np.kron(sy,I1),sy))
    H += Jz2[i]    *         np.kron(np.kron(sz,I1),sz)
    
   
    H += hz[i]   * np.kron(np.kron(sz,I1),I1) +  hx[i]   * np.kron(np.kron(sx,I1),I1)
    H += hz[i+1] * np.kron(np.kron(I1,sz),I1) +  hx[i+1] * np.kron(np.kron(I1,sx),I1)
    H += hz[i+2] * np.kron(np.kron(I1,I1),sz) +  hx[i+2] * np.kron(np.kron(I1,I1),sx)

    H_bond.append(np.reshape(H,(d,d,d,d,d,d)))
    expH_bond.append(np.reshape(expm(-delta*H),(d,d,d,d,d,d)))
        
    
    return expH_bond,H_bond


###############################################
#       Local Hamiltonian infinite algorithm
###############################################
def mixing_hamiltinian_expH_bond(Jxy1,Jz1,Jxy2,Jz2,hz,hx,L,delta):
    """ Returns bond hamiltonian and bond time-evolution"""
    sx = np.array([[0., 1.],[1., 0.]])
    sy = np.array([[0.,-1j],[1j, 0.]])
    sz = np.array([[1., 0.],[0.,-1.]])
    d = 4
    
    I1 = np.eye(2)
     

    expH_bond = []
    H_bond = []
    
    for i in range(0,L-2,2):
    
        # first neigbour
        H  =  Jxy1[i]  * np.real(np.kron(np.kron(np.kron(sx,sx),I1),I1) + np.kron(np.kron(np.kron(sy,sy),I1),I1))
        H +=   Jz1[i]  *         np.kron(np.kron(np.kron(sz,sz),I1),I1)
        
        H += Jxy1[i+1] * np.real(np.kron(np.kron(np.kron(I1,sx),sx),I1) + np.kron(np.kron(np.kron(I1,sy),sy),I1))
        H +=  Jz1[i+1] *         np.kron(np.kron(np.kron(I1,sz),sz),I1)
        
        # secend neigbour
        H +=  Jxy2[i]  * np.real(np.kron(np.kron(np.kron(sx,I1),sx),I1) + np.kron(np.kron(np.kron(sy,I1),sy),I1))
        H +=   Jz2[i]  * np.kron(np.kron(np.kron(sz,I1),sz),I1)
        
        H += Jxy2[i+1] * np.real(np.kron(np.kron(np.kron(I1,sx),I1),sx) + np.kron(np.kron(np.kron(I1,sy),I1),sy))
        H +=  Jz2[i+1] * np.kron(np.kron(np.kron(I1,sz),I1),sz)
        
        # onsite-term
        H +=     hz[i] * np.kron(np.kron(np.kron(sz,I1),I1),I1)
        H +=     hx[i] * np.kron(np.kron(np.kron(sx,I1),I1),I1)
        
        H +=   hz[i+1] * np.kron(np.kron(np.kron(I1,sz),I1),I1)
        H +=   hx[i+1] * np.kron(np.kron(np.kron(I1,sx),I1),I1)
        
        H_bond.append(np.reshape(H,(d,d,d,d)))
        expH_bond.append(np.reshape(expm(-delta*H),(d,d,d,d)))
    
    
    i = L-2
    # first neigbour
    H  = Jxy1[i] * np.real( np.kron(np.kron(np.kron(sx,sx),I1),I1) + np.kron(np.kron(np.kron(sy,sy),I1),I1))
    H +=  Jz1[i] * np.kron(np.kron(np.kron(sz,sz),I1),I1)
    
    # onsite-term
    H +=   hz[i] * np.kron(np.kron(np.kron(sz,I1),I1),I1)
    H +=   hx[i] * np.kron(np.kron(np.kron(sx,I1),I1),I1)
    H += hz[i+1] * np.kron(np.kron(np.kron(I1,sz),I1),I1)
    H += hx[i+1] * np.kron(np.kron(np.kron(I1,sx),I1),I1)
    
        

    return expH_bond,H_bond

###############################################
#       Local Hamiltonian infinite algorithm
###############################################
def infinite_hamiltinian_expH_bond(Jxy,Jz,hz,hx,L,delta,infinite=False):
    """ Returns bond hamiltonian and bond time-evolution"""
    sx = np.array([[0., 1.],[1., 0.]])
    sy = np.array([[0.,-1j],[1j, 0.]])
    sz = np.array([[1., 0.],[0.,-1.]])
    d = 2

    expH_bond = []
    H_bond = []
    
    for i in range(L):
    	H = Jxy[i]*(np.kron(sx,sx) + np.kron(sx,sx)) + Jz[i]*np.kron(sz,sz)
    	H = H + hz[i]*np.kron(sz,np.eye(2))
    	H = H + hx[i]*np.kron(sx,np.eye(2))
    	H_bond.append(np.reshape(H,(d,d,d,d)))
    	expH_bond.append(np.reshape(expm(-delta*H),(d,d,d,d)))

    return expH_bond,H_bond
    
    
  
###############################################
#       Local Hamiltonian nearest neighbour
###############################################
def pott_expH_bond(J,f,L,delta):
    """ Returns bond hamiltonian and bond time-evolution"""
    d = 3
    w = np.exp(-1j*2*np.pi/d)
    sig = np.array([[1., 0., 0.],[0., w , 0.],[0., 0., w**2]])
    tau = np.array([[0., 0., 1.],[1., 0., 0.],[0., 1.,   0.]])
    
    I3 = np.eye(3)

    expH_bond = []
    H_bond = []
   
    for i in range(L-2):
        H  =  J[i] * np.kron(sig.conj(), sig        )
        H +=  J[i] * np.kron(sig       , sig.conj() )
        H +=  f[i] * np.kron(tau       , I3         )
        H +=  f[i] * np.kron(tau.T     , I3         )
        H +=  f[i] * np.kron(I3        , tau        )
        H +=  f[i] * np.kron(I3        , tau.T      )
        H_bond.append(np.reshape(H,(d,d,d,d)))
        expH_bond.append(np.reshape(expm(-delta*H),(d,d,d,d)))

    i  =  L-2
    H  = J[i] * np.kron(sig.conj(),        sig )
    H += J[i] * np.kron(sig       , sig.conj() )
    
    H += f[i] * np.kron(tau , I3  ) + f[i+1] * np.kron(tau , I3  )
    H += f[i] * np.kron(I3  , tau ) + f[i+1] * np.kron(I3  , tau )

    H_bond.append(np.reshape(H,(d,d,d,d)))
    expH_bond.append(np.reshape(expm(-delta*H),(d,d,d,d)))
        
    return expH_bond,H_bond
########################################
#       test AB
########################################
def testAB(l, G):
    A, B = lG_AB(l, G)
    L    = len(G)
    
    for pos in range(L):
        testA = np.tensordot( np.conj(A[pos]), A[pos], axes=([0,1], [0,1]))
        testB = np.tensordot( B[pos], np.conj(B[pos]), axes=([2,0], [2,0]))
    
        if np.sum(np.abs( testA - np.eye(testA.shape[0]) )) / np.sum(np.abs(A[pos])) > 1.0e-10 : print(' left normalize')
        if np.sum(np.abs( testB - np.eye(testB.shape[0]) )) / np.sum(np.abs(B[pos])) > 1.0e-10 : print(' right normalize')


########################################
#    Left-Normalization for Vidal MPS
########################################
def normalize_left(l, G):
    A, B = lG_AB(l, G)
    AB = B
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
           linvU    = np.tensordot ( np.diag(1.0/(l[pos]+1e-20)), U              , axes=(1,1) )
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
    AB = A
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
           Vlinv    = np.tensordot ( V     , np.diag(1.0/(l[pos+1]+1e-20)), axes=(2,0) )
           G[pos]   = Vlinv

        if pos > 0 :  l[pos] = S
        else       :  norm = (np.tensordot( U, np.diag(S), axes=(1,0))[0,0]).real
    return l, G, norm
           

########################################
#    test Left Or Right Normalie
########################################

def testLeftOrRightNormalie(l, G):
    L = len(G)
    A, B = lG_AB(l, G)
    testA = []
    testB = []
    for pos in range(L):
        testA.append(abs(np.tensordot( np.conj(A[pos]) ,  A[pos]           , axes=([0,1,2],[0,1,2]) )))
        testB.append(abs(np.tensordot( B[pos]          ,  np.conj(B[pos])  , axes=([0,1,2],[0,1,2]) )))
    if all(abs(1.0-i)<1e-6 for i in testA):
        #print('left  normalize')
        return 1
    if all(abs(1.0-i)<1e-6 for i in testB):
        #print('right normalize')
        return 1
        
########################################
#           MPS overlap
########################################

def overlap(l1, G1, l2, G2):
    L = len(G1)
    A1, B1 = lG_AB(l1, G2)
    A2, B2 = lG_AB(l1, G2)
            
    #if testLeftOrRightNormalie()
    tem = 1
    for pos in range(L):
        tem *= np.abs(np.tensordot( np.conj(A1[pos]),   A2[pos] , axes=([0,1,2],[0,1,2])))
    return tem


def overlap_v2(l1, G1, l2, G2):
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
    


########################################
#          tensor product of two op
########################################
def op_tensor(op1,op2):
    top = np.tensordot( op1[:,np.newaxis], op2[np.newaxis,:], axes=(1,0) )
    top = np.transpose( top              , (0,2,1,3)                     )
    return top

def op_tensor_v2(op1,op2):
    d = 2
    return np.reshape(np.kron(op1,op2),(d,d,d,d))

########################################
#           one-site expection
########################################
def expc_one_site(l,G,op,pos):
    ''' expetion of one site operator at : pos '''
    Psi = np.tensordot( np.diag(l[pos]),  G[pos]            , axes=(1,1)         )
    Psi = np.tensordot( Psi            ,  np.diag(l[pos+1]) , axes=(2,0)         )
    Psi = np.tensordot( np.conj(Psi)   ,  Psi               , axes=([0,2],[0,2]) )
    res = np.tensordot( Psi            ,  op                , axes=([0,1],[0,1]) )
    return res
    
    
########################################
#    two-sites expection (nearest-neighbor)
########################################
def expc_two_site(l,G,op,pos):
    ''' expetion of two site operator at : pos,pos+1 '''
    Psi = np.tensordot( np.diag(l[pos]),  G[pos]            , axes=(1,1)                 )
    Psi = np.tensordot( Psi            ,  np.diag(l[pos+1]) , axes=(2,0)                 )
    Psi = np.tensordot( Psi            ,  G[pos+1]          , axes=(2,1)                 )
    Psi = np.tensordot( Psi            ,  np.diag(l[pos+2]) , axes=(3,0)                 )
    Psi = np.tensordot( np.conj(Psi)   ,  Psi               , axes=([0,3],[0,3])         )
    res = np.tensordot( Psi            ,  op[pos]           , axes=([0,1,2,3],[0,1,2,3]) )
    return res


    
########################################
#       Enatanglement entropy
########################################
def entanglement_entropy(l):
    " Returns the half chain entanglement entropy "
    S=[]
    L = len(l)
    for i_bond in range(L):
        x=l[i_bond][l[i_bond]>10**(-20)]**2
        S.append(-np.inner(np.log(x),x))
    return S
    

########################################
#       Enatanglement entropy
########################################
def loschmidt_echo(l0,G0,l,G):
    " loschmidt_echo  "
    
    E=[]
    L = len(G)
    for pos in range(L-1):
        sG0 = np.tensordot( np.diag(l0[pos]) , G0[pos]   , axes=(1,1))
        C   = np.tensordot( sG0              , np.eye(2) , axes=(1,0))
        sG  = np.tensordot( np.diag(l[pos])  , G[pos]    , axes=(1,1))
        sG  = sG.conj()
        E.append(np.squeeze(np.tensordot(sG,C,axes=([0,1,2],[0,2,1]))).item())
    return E
    
def loschmidt_echo_v2(l0,G0,l,G):
    " loschmidt_echo  "
    return overlap(l0,G0,l,G)
    
    
    
###########################################################
#     First order Suzuki-Trotter decomposition time evolve
##########################################################
def Suzuki_Trotter_1ftOrder_two_sites(l, G, expH, infinite=False):
    """ Perform the real or imaginary time evolution of the MPS """
    L = len(G)
    
    for k in [0,1]:
        for i_bond in range(k, L-1, 2):
            i0 = i_bond
            if infinite:
                i1 = np.mod(i_bond+1,L)
                i2 = np.mod(i_bond+2,L)
            else:
                i1 = i_bond+1
                i2 = i_bond+2
            
            d     = G[i0].shape[0]
            chi0  = G[i0].shape[1]
            chi1  = G[i0].shape[2]
            chi2  = G[i1].shape[2]
            
            l0inv = 1.0/(l[i0]+1e-20)
            l2inv = 1.0/(l[i2]+1e-20)
            
            # step A
            Psi = np.tensordot( np.diag(l[i0]), G[i0]          , axes=(1,1) )
            Psi = np.tensordot( Psi           , np.diag(l[i1]) , axes=(2,0) )
            Psi = np.tensordot( Psi           , G[i1]          , axes=(2,1) )
            Psi = np.tensordot( Psi           , np.diag(l[i2]) , axes=(3,0) )
            
            # step B
            Phi = np.tensordot( Psi, expH[i0], axes=([1,2],[0,1]) )
            Phi = np.transpose( Phi, (2,0,3,1)                    )
            Phi = np.reshape  ( Phi, (d*chi0,d*chi2)              )
            
            # step C
            U, S, V = svd( Phi, full_matrices=0 )
            disc = np.sum( S[chi1:]**2 ) / np.sum( S**2 )
            
            
            S = S[0:chi1] / np.sqrt( np.sum(S[0:chi1]**2) )                    # new Lambda1
            l[i1] = S
            
            U = np.reshape  ( U[:,0:chi1], (d,chi0,chi1)              )
            U = np.tensordot( U          , np.diag(l0inv), axes=(1,0) )
            U = np.transpose( U          , (0,2,1)                    )        # new Gamma0
            G[i0] = U
            
            V = np.reshape  ( V[0:chi1,:], (chi1,d,chi2)              )
            V = np.transpose( V          , (1,0,2)                    )
            V = np.tensordot( V          , np.diag(l2inv), axes=(2,0) )        # new Gamma1
            G[i1] = V
            

###########################################################
#     First order Suzuki-Trotter decomposition time evolve
##########################################################
def Suzuki_Trotter_1ftOrder_three_sites(l, G, expH,infinite=False):
    """ Perform the real or imaginary time evolution of the MPS """
    L = len(G)
    if np.mod(L,3)!=0:
        sys.exit('ERR: # np.mod(L,3)!=0 ')
    
    for k in [0,1]:
        for i_bond in range(k, L-2, 3):
            i0 = i_bond
            i1 = i_bond+1
            i2 = i_bond+2
            i3 = i_bond+3
            
            d     = G[i0].shape[0]
            chi0  = G[i0].shape[1]
            chi1  = G[i0].shape[2]
            chi2  = G[i1].shape[2]
            chi3  = G[i2].shape[2]
            
            l0inv = 1.0/(l[i0]+1e-20)
            #l1inv = 1.0/(l[i1]+1e-20)
            l3inv = 1.0/(l[i3]+1e-20)
            
            # step A
            Psi = np.tensordot( np.diag(l[i0]), G[i0]          , axes=(1,1) )
            Psi = np.tensordot( Psi           , np.diag(l[i1]) , axes=(2,0) )
            Psi = np.tensordot( Psi           , G[i1]          , axes=(2,1) )
            Psi = np.tensordot( Psi           , np.diag(l[i2]) , axes=(3,0) )
            Psi = np.tensordot( Psi           , G[i2]          , axes=(3,1) )
            Psi = np.tensordot( Psi           , np.diag(l[i3]) , axes=(4,0) )
            
            # step B
            Phi = np.tensordot( Psi, expH[i0] , axes=([1,2,3],[0,1,2]) )
            Phi = np.transpose( Phi, (2,0,3,4,1)                       )
            Phi = np.reshape  ( Phi, (d*chi0,d*d*chi3)                 )
            
            # step C first SVD
            U, S, V = svd( Phi, full_matrices=0 )
            disc1 = np.sum( S[chi1:]**2 ) / np.sum( S**2 )
            
            
            S = S[0:chi1] / np.sqrt( np.sum(S[0:chi1]**2) )                    # new Lambda1
            l[i1] = S
            l1inv = 1.0/(l[i1]+1e-20)
            
            U = np.reshape  ( U[:,0:chi1]   , (d,chi0,chi1)               )
            U = np.tensordot( np.diag(l0inv), U            ,   axes=(1,1) )
            U = np.transpose( U             , (1,0,2)                     )
            G[i0] = U                                                          # new Gamma0
            
            Phi = np.reshape  ( V[0:chi1,:]   , (chi1,d,d,chi3)             )
            Phi = np.tensordot( np.diag(l1inv), Phi            , axes=(1,0) )
            Phi = np.transpose( Phi           , (1,0,2,3)                   )  
            Phi = np.reshape  ( Phi           , (d*chi1,d*chi3)             )
            
            
            # step D second SVD
            U, S, V = svd( Phi, full_matrices=0 )
            disc2 = np.sum( S[chi2:]**2 ) / np.sum( S**2  )
            
            S = S[0:chi2] / np.sqrt( np.sum(S[0:chi2]**2) )
            l[i2] = S                                                          # new Lambda2
            
            
            U = np.reshape  ( U[:,0:chi2]   , (d,chi1,chi2)               )
            U = np.tensordot( np.diag(l1inv), U            ,   axes=(1,1) )
            U = np.transpose( U             , (1,0,2)                     )
            G[i1] = U                                                          # new Gamma1
            
            
            V = np.reshape  ( V[0:chi2,:], (chi2,d,chi3)              )
            V = np.transpose( V          , (1,0,2)                    )
            V = np.tensordot( V          , np.diag(l3inv), axes=(2,0) )
            G[i2] = V                                                          # new Gamma2
            
            
            
            
        


     
     
