from scipy.linalg import expm 
import pylab as pl
import numpy as np
from svd_robust import svd

""" Conventions:
B[i,a,b] has axes (physical, left virtual, right virtual),
W[a,b,i,j] has axes (virtual left, virtual right, physical out, physical in)
S[i] are schmidt values between sites (i, i+1),
H_bond[i] is the bond hamiltonian between (i,i+1) with (only physical)
axes (out left, out right, in left, in right)"""

"""
def init_af_mps(L):
    ''' Returns FM Ising MPS
'''
    d = 2
    B = []
    s = []
    for i in range(L):
        B.append(np.zeros([d,1,1]))
        B[-1][np.mod(i,d),0,0]=1
        #B[-1][i%d,0,0]=1
        s.append(np.ones([1])) # set trivial initial weights
    s.append(np.ones([1]))     # set trivial initial weights
    return B,s
"""

def init_af_mps(L):
    chiMax = 20
    d = 2
    B = []
    s = []
    chi = [0] * (L+1)

    for pos in range(0,L//2+1):
        chi[pos] = np.minimum( chiMax , d**pos     )
    for pos in range(L//2+1,L+1):
        chi[pos] = np.minimum( chiMax , d**(L-pos) )

    for pos in range(L):
        chi0 = int( chi[ pos    ] )
        chi1 = int( chi[ pos+1  ] )
        B.append( np.zeros( (d,chi0,chi1), np.complex128 ) )
        s.append( np.zeros( (  chi0     )                ) )

    chi0 = int( chi[ L    ] )
    s.append(np.zeros( (  chi0     )                ))
    ########################################
    #       Initial state
    ########################################
    for pos in range(L):
        B[pos][pos%2,0,0] = 1.0 # example: Neel state
        s[pos][   :     ] = 1.0
    s[L][   :     ] = 1.0
    return B, s




def init_heisenberg_U_bond(J,D,h,L,delta):
    """ Returns bond hamiltonian and bond time-evolution"""
    sx = np.array([[0., 1.],[1., 0.]])/2.
    sy = np.array([[0.,-1j],[1j, 0.]])/2.
    sz = np.array([[1., 0.],[0.,-1.]])/2.
    d = 2
    
    U_bond = []
    H_bond = []
    for i in range(L-2):
        H = J[i]*np.real((np.kron(sx,sx) + np.kron(sy,sy))) + D[i]*np.kron(sz,sz) 
        H = H + h[i]*np.kron(sz,np.eye(2))
        H_bond.append(np.reshape(H,(d,d,d,d)))
        U_bond.append(np.reshape(expm(-delta*H),(d,d,d,d)))
    
    i = L-2
    H = J[i]*np.real((np.kron(sx,sx) + np.kron(sy,sy))) + D[i]*np.kron(sz,sz) 
    H = H + h[i]*np.kron(sz,np.eye(2)) + h[i+1]*np.kron(np.eye(2),sz)
    
    H_bond.append(np.reshape(H,(d,d,d,d)))
    U_bond.append(np.reshape(expm(-delta*H),(d,d,d,d)))
    return U_bond,H_bond

def entanglement_entropy(s):
    " Returns the half chain entanglement entropy "
    S=[]
    L = len(s)
    for i_bond in range(L):
        x=s[i_bond][s[i_bond]>10**(-20)]**2
        S.append(-np.inner(np.log(x),x))
    return(S)
	
 
def site_expectation(B,s,O_list):
    " Expectation value for a site operator "
    E=[]
    L = len(B)
    for isite in range(L-1):
        sB = np.tensordot(np.diag(s[isite]),B[isite],axes=(1,1))
        C = np.tensordot(sB,O_list[isite],axes=(1,0))
        sB=sB.conj()
        E.append(np.squeeze(np.tensordot(sB,C,axes=([0,1,2],[0,2,1]))).item())
    return(E)
    
    
def bond_expectation(B,s,H_bond):
    " Expectation value for a bond operator "
    E=[]
    L = len(B)
    for i_bond in range(L-1):
        BB = np.tensordot(B[i_bond],B[i_bond+1],axes=(2,1))
        sBB = np.tensordot(np.diag(s[i_bond]),BB,axes=(1,1))
        C = np.tensordot(sBB,H_bond[i_bond],axes=([1,2],[2,3]))
        sBB=np.conj(sBB)
        E.append(np.squeeze(np.tensordot(sBB,C,axes=([0,3,1,2],[0,1,2,3]))).item()) 
    return E


    
    
def loschmidt_echo(B0,s0,B,s):
    " loschmidt_echo  "
    E=[]
    L = len(B)
    for isite in range(L-1):
        sB0 = np.tensordot(np.diag(s0[isite]),B0[isite],axes=(1,1))
        C = np.tensordot(sB0,np.eye(2),axes=(1,0))
        sB = np.tensordot(np.diag(s[isite]),B[isite],axes=(1,1))
        sB=sB.conj()
        E.append(np.squeeze(np.tensordot(sB,C,axes=([0,1,2],[0,2,1]))).item())
    return(E)

def sweep(B,s,U_bond,chimax,eps=10e-8):
    """ Perform the imaginary time evolution of the MPS """
    L = len(B)
    d = B[0].shape[0]
    for k in [0,1]:
        for i_bond in range(k, L-1,2):
            ia = i_bond
            ib = i_bond+1
            ic = i_bond+2
            chia = B[ia].shape[1]
            chic = B[ib].shape[2]

            
            # Construct theta matrix and time evolution #
            theta = np.tensordot(B[ia],B[ib],axes=(2,1))                      # i  a  j  b
            theta = np.tensordot(U_bond[i_bond],theta,axes=([2,3],[0,2]))     # ip jp a  b
            theta = np.tensordot(np.diag(s[ia]),theta,axes=([1,2]))           # a  ip jp b
            theta = np.reshape(np.transpose(theta,(1,0,2,3)),(d*chia,d*chic)) # ip a  jp b

            # Schmidt deomposition #
            X, Y, Z = svd(theta,full_matrices=0)
            chi2 = np.min([np.sum(Y>eps), chimax])

            piv = np.zeros(len(Y), np.bool)
            piv[(np.argsort(Y)[::-1])[:chi2]] = True

            Y = Y[piv]; invsq = np.sqrt(sum(Y**2))
            X = X[:,piv]
            Z = Z[piv,:]

            # Obtain the new values for B and s #
            s[ib] = Y/invsq

            X = np.reshape(X,(d,chia,chi2))
            #back to right canonical form
            X = np.transpose(np.tensordot(np.diag(s[ia]**(-1)),X,axes=(1,1)),(1,0,2))
            B[ia] = np.tensordot(X, np.diag(s[ib]),axes=(2,0))
            B[ib] = np.transpose(np.reshape(Z,(chi2,d,chic)),(1,0,2))
            
            """
            
            # construct theta
            theta = np.tensordot(np.diag(s[ia]),B[ia] , axes=(1,1))
            theta = np.tensordot(theta, np.diag(s[ib]), axes=(2,0))
            theta = np.tensordot(theta, B[ib]         , axes=(2,1))
            theta = np.tensordot(theta, np.diag(s[ic]), axes=(3,0))
            
            # Apply U
            theta = np.tensordot(theta,U_bond[i_bond],axes=([1,2],[0,1]))
            theta = np.transpose(theta,(2,0,3,1))
            theta = np.reshape(theta,(d*chia,d*chic))
            
            # SVD
            X, Y, Z = svd(theta,full_matrices=0)
            chi2 = np.min([np.sum(Y>eps), chimax])
            
            Y = Y[0:chi2]; invsq = np.sqrt(sum(Y**2))
            X = X[:,0:chi2]
            Z = Z[0:chi2,:]
             
            # Obtain the new values for B and s #
            s[ib] = Y/invsq 

            X = np.reshape(X,(d,chia,chi2))
            X = np.tensordot( X, np.diag(s[ia]**(-1)),axes=(1,0))
            B[ia] = np.transpose(X,(0,2,1))
            
            
            Z = np.reshape(Z,(chi2,d,chic))
            Z = np.transpose(Z,(1,0,2))
            B[ib] = np.tensordot( Z, np.diag(s[ic]**(-1)),axes=(2,0))
            """
    

def clean_quench():
    ######## Define the simulation parameter ######################
    chi = 50

    ######## Define the model parameter ###########################
    L = 40
    J = (L-1)*[1.]
    D = (L-1)*[1.]
    h = L*[0.]
    delta = 0.1
    T = 10
    N = int(T/delta)
    sz = np.array([[1.,0.],[0.,-1.]])/2.
    
    ######### Entanglement growth ##########
    B , s = init_af_mps(L)
    B0 =  B 
    s0 = s
    
    U_bond,H_bond = init_heisenberg_U_bond(J,D,h,L,delta*1j)
    S = []
    m = []
    LE = []
    bond = []
    t = []
    for i in range(N):
        
        sweep(B,s,U_bond,chi)
        t.append((i+1)*delta)
        S.append(entanglement_entropy(s)[L//2])
        #m.append(np.abs(site_expectation(B,s,L*[sz])[L//2]))
        m.append(np.sum(site_expectation(B,s,L*[sz])))
        bond.append(np.sum(bond_expectation(B,s,H_bond)))
        LE.append(np.abs(np.sum(loschmidt_echo(B0,s0, B,s)))**2)
        
   

    ######## Plot ##########################
    #pl.plot(t,S)
    #pl.ylim([0,np.max(S)])
    #pl.ylabel('$S$')
    #pl.xlabel('$t$')
    
    pl.figure()
    pl.plot(t,m)
    pl.ylim([np.min(m)-0.1,np.max(m)+0.1])
    pl.ylabel('$m$')
    pl.xlabel('$t$')
    pl.show()
    
    #pl.figure()
    #pl.plot(t,LE)
    #pl.ylim([0,np.max(LE)])
    #pl.ylabel('$LE$')
    #pl.xlabel('$t$')
    #pl.show()



def disordered_quench():    
    ######## Define the simulation parameter ######################
    chi = 20

    ######## Define the model parameter ###########################
    L = 20
    J = (L-1)*[1.]
    D = (L-1)*[1.]
    W = 8.
    N_disorder = 20
    delta = 0.1
    T = 20
    N = int(T/delta)
    sz = np.array([[1.,0.],[0.,-1.]])/2.

    ######### Entanglement growth ##########
    S_list = []
    m_list = []

    for disorder in range(N_disorder):
        print(disorder)
        h = 2*W*(0.5-np.random.rand(L))
        U_bond,H_bond = init_heisenberg_U_bond(J,D,h,L,delta*1j)
    
        S = []
        m = []
        t = []
        B,s = init_af_mps(L)
    
        for i in range(N):
            sweep(B,s,U_bond,chi)
            t.append((i+1)*delta)
            S.append(entanglement_entropy(s)[L//2])
            m.append(np.abs(site_expectation(B,s,L*[sz])[L//2]))
            
        S_list.append(S)
        m_list.append(m)

    ######## Plot ##########################
    pl.plot(t,np.mean(S_list,axis=0))
    pl.ylabel('$S$')
    pl.xlabel('$t$')
    pl.ylim([0,np.max(S_list)])
    pl.figure()
    pl.plot(t,np.mean(m_list,axis=0))
    pl.ylim([0,np.max(m_list)])
    pl.ylabel('$m$')
    pl.xlabel('$t$')
    pl.show()
   

clean_quench()
#disordered_quench()
