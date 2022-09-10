import numpy as np
from scipy import linalg
from svd_robust import svd
from scipy.linalg import expm, sinm, cosm
import sys,os

##############################
#        inverse getting rid of zero devision
###############################
def inverse(ll):
    a = [ll>1e-10]
    b = 1.0/ll[a[0]]
    c = np.zeros(len(ll)-len(b))
    return np.concatenate((b,c))

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
    return l, G, chi

###############################################
#       Local Hamiltonian nearest neighbour
###############################################
def hamiltinian(Jxy,Jz,hz,hx):
    """ Returns bond hamiltonian and bond time-evolution"""
    sx = np.array([[0., 1.],[1., 0.]])
    sy = np.array([[0.,-1j],[1j, 0.]])
    sz = np.array([[1., 0.],[0.,-1.]])
    d  = 2
    I1 = np.eye(2)
    H  = Jxy * np.real(np.kron(sx,sx) + np.kron(sy,sy))
    H += Jz  *         np.kron(sz,sz)
    H += hz  *         np.kron(sz,I1)
    H += hx  *         np.kron(sx,I1)
    H += hz  *         np.kron(I1,sz)
    H += hx  *         np.kron(I1,sx)

    return np.reshape(H, (d,d,d,d))   
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

    
def overlap_inf(l1, G1, l2, G2):
    L = len(G1)
    A1, B1 = lG_AB(l1, G1)
    A2, B2 = lG_AB(l2, G2)

        
    #pos = 0
    #OL_now =  np.tensordot (np.conj(A1[pos]),  A2[pos]          , axes=([0,1],[0,1]) )
    for pos in range(L):
        pos    = np.mod(pos,L)
        OL_now =  np.tensordot(np.conj(A1[pos]),  A2[pos]          , axes=([0,1],[0,1]) )
        pos    = np.mod(pos+1,L) # ib==pos
        temp   = np.tensordot(OL_now           ,  np.conj(A1[pos]) , axes=(0,1)         )
        temp   = np.tensordot(temp             ,  A2[pos]          , axes=([1,0],[0,1]) )
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
def expvalue1(l, G, op, pos):
    ''' expetion of one site operator at : pos '''
    Psi = np.tensordot( np.diag(l[pos]),  G[pos]            , axes=(1,1)         )
    Psi = np.tensordot( Psi            ,  np.diag(l[pos+1]) , axes=(2,0)         )
    Psi = np.tensordot( np.conj(Psi)   ,  Psi               , axes=([0,2],[0,2]) )
    print(ia,ib,Psi.shape,op[ib].shape)
    res = np.tensordot( Psi            ,  op                , axes=([0,1],[0,1]) )
    return res
    
def expvalue1_ibc(l, G, op, pos):
    " Expectation value for a site operator "
    L = len(G)
    ia  = np.mod(pos,L)
    ib  = np.mod(pos+1,L) # ib==pos
    Psi = np.tensordot( np.diag(l[ib])   ,  G[ib]          , axes=(1,1)         )
    Psi = np.tensordot( Psi              ,  np.diag(l[ia]) , axes=(2,0)         )
    Psi = np.tensordot( np.conj(Psi)     ,  Psi            , axes=([0,2],[0,2]) )
    res = np.tensordot( Psi              ,  op             , axes=([0,1],[0,1]) )
    return res

def expvalue1_inf(l0, G0, l, G, op, pos):
    " Expectation value for a site operator "
    L = len(G)
    ia  = np.mod(pos,L)
    ib  = np.mod(pos+1,L) # ib==pos
    Psi0 = np.tensordot( np.diag(l0[ib])  ,  G0[ib]         , axes=(1,1)         )
    Psi0 = np.tensordot( Psi0             ,  np.diag(l0[ia]), axes=(2,0)         )
    Psi  = np.tensordot( np.diag(l[ib])   ,  G[ib]          , axes=(1,1)         )
    Psi  = np.tensordot( Psi              ,  np.diag(l[ia]) , axes=(2,0)         )
    
    Psi = np.tensordot( np.conj(Psi0)     ,  Psi            , axes=([0,2],[0,2]) )
    res = np.tensordot( Psi               ,  op             , axes=([0,1],[0,1]) )
    return res
########################################
#    two-sites expection (nearest-neighbor)
########################################
def expvalue2(l, G, op, pos):
    ''' expetion of two site operator at : pos,pos+1 '''
    Psi = np.tensordot( np.diag(l[pos]),  G[pos]            , axes=(1,1)                 )
    Psi = np.tensordot( Psi            ,  np.diag(l[pos+1]) , axes=(2,0)                 )
    Psi = np.tensordot( Psi            ,  G[pos+1]          , axes=(2,1)                 )
    Psi = np.tensordot( Psi            ,  np.diag(l[pos+2]) , axes=(3,0)                 )
    Psi = np.tensordot( np.conj(Psi)   ,  Psi               , axes=([0,3],[0,3])         )
    res = np.tensordot( Psi            ,  op                , axes=([0,1,2,3],[0,1,2,3]) )
    return res

def expvalue2_ibc(l, G, op, pos):
    " Expectation value for a bound operators "
    L = len(G)
    ia  = np.mod(pos,L)
    ib  = np.mod(pos+1,L)
    Psi = np.tensordot( np.diag(l[ib])  ,  G[ia]          , axes=(1,1)                 )
    Psi = np.tensordot( Psi             ,  np.diag(l[ia]) , axes=(2,0)                 )
    Psi = np.tensordot( Psi             ,  G[ib]          , axes=(2,1)                 )
    Psi = np.tensordot( Psi             ,  np.diag(l[ib]) , axes=(3,0)                 )
    Psi = np.tensordot( np.conj(Psi)    ,  Psi            , axes=([0,3],[0,3])         )
    res = np.tensordot( Psi             ,  op             , axes=([0,1,2,3],[0,1,2,3]) )
        
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

##############################
#        time_evolution
###############################
def time_evolution(l,G,H,BC,delta,numsteps):
    L = len(G)
    d = G[0].shape[0]
    
    if BC=='OBC' :
       #  bond 1,2,... : between sites 0-1, 1-2,...
       evenbonds = np.arange(2,L,2)
       oddbonds  = np.arange(1,L,2)
    else: 
       evenbonds = np.arange(0,L,2)
       oddbonds  = np.arange(1,L,2)

    H      = np.reshape(H,     (d*d,d*d) )
    expH1  =       expm(       delta * H )
    expH2  =       expm( 0.5 * delta * H )
    expH1  = np.reshape(expH1, (d,d,d,d) )
    expH2  = np.reshape(expH2, (d,d,d,d) ) 

    disc = 0.0; temp=0.0
    l, G, tem = evolve_bonds(l,G,expH2,evenbonds ); disc += temp

    #for num in range(numsteps-1):
    #    l, G, tem = evolve_bonds(l,G,expH1,oddbonds ); disc += temp
    #    l, G, tem = evolve_bonds(l,G,expH1,evenbonds); disc += temp

    l, G, tem = evolve_bonds(l,G,expH1,oddbonds  ); disc += temp
    l, G, tem = evolve_bonds(l,G,expH2,evenbonds ); disc += temp
    
    return l, G, disc



##############################
#        evolve_bonds
###############################
def evolve_bonds(l,G,expH,bonds):
    L    = len(G)
    disc = 0.0
    for pos in bonds:
        posM  = np.mod(pos-1,L)
        posP  = np.mod(pos+1,L)
        G[posM], l[pos], G[pos], temp = evolve_single_bond(l[posM], G[posM], l[pos], G[pos], l[posP],expH)
        disc += temp
        
    return l, G, disc


###########################################################
#     evolve_single_bond
##########################################################
def evolve_single_bond(l0, G0, l1, G1, l2, expH):
    """ Perform the real or imaginary time evolution of the MPS """
   
    d     = G0.shape[0]
    chi0  = G0.shape[1]
    chi1  = G0.shape[2]
    chi2  = G1.shape[2]
    
    l0inv = inverse(l0)# or 1.0/(l0+1e-20)
    l2inv = inverse(l2)# or 1.0/(l2+1e-20)
    
    # step A
    Psi = np.tensordot( np.diag(l0), G0          , axes=(1,1) )
    Psi = np.tensordot( Psi        , np.diag(l1) , axes=(2,0) )
    Psi = np.tensordot( Psi        , G1          , axes=(2,1) )
    Psi = np.tensordot( Psi        , np.diag(l2) , axes=(3,0) )
    
    # step B
    Phi = np.tensordot( Psi, expH, axes=([1,2],[0,1]) )
    Phi = np.transpose( Phi, (2,0,3,1)                )
    Phi = np.reshape  ( Phi, (d*chi0,d*chi2)          )
            
    # step C
    U, S, V = svd( Phi, full_matrices=0 )
    disc = np.sum( S[chi1:]**2 ) / np.sum( S**2 )
              
    l1  = S[0:chi1] / np.sqrt( np.sum(S[0:chi1]**2) )                   # new Lambda1
          
    U  = np.reshape  ( U[:,0:chi1], (d,chi0,chi1)              )
    U  = np.tensordot( U          , np.diag(l0inv), axes=(1,0) )
    G0 = np.transpose( U          , (0,2,1)                    )        # new Gamma0
                    
    V  = np.reshape  ( V[0:chi1,:], (chi1,d,chi2)              )
    V  = np.transpose( V          , (1,0,2)                    )
    G1 = np.tensordot( V          , np.diag(l2inv), axes=(2,0) )        # new Gamma1
          
    return G0, l1, G1, disc  

