import scipy as sp
import scipy.sparse as sparse
import numpy as np
import math
from numba import njit

#@njit
def SvN(state,i,N):
    
    #%%
    ###################################
    #  Von Numann entropy with
    #  Singula value decomposition (SVD)
    ####################################
    
    
    gs = np.reshape(state,(int(2**i),int(2**(N-i))))
    
    s = sparse.linalg.svds(gs,
                              min(gs.shape)-2,
                              which='LM',
                              return_singular_vectors=False)
           
    p=abs(s)**2
    p[p <= 0] = 1e-20
    EE_svd = -np.nansum(p*np.log(p))
    
    svn  = -np.sum(p*np.log(p))
    if math.isnan(svn):
        return  svn
    
    return  EE_svd
    '''
    #%%
    ###################################
    #  Von Numann entropy with
    #  Reduced Density Matrix (RDM)
    ####################################
    da = 2**(N-i)
    db = 2**i
    rho = np.reshape(state,(da,db)).conj().T @ np.reshape(state,(da,db))
    eigs = sp.linalg.eigh(rho,eigvals_only=True)
    del(rho) # clear memory
    
    eigs[eigs <= 0] = 1e-20
    svn  = -np.sum(np.real(eigs*np.log(eigs)))
    if math.isnan(svn):
        EE_rdm = svn
    else:
        EE_rdm = svn
    
    return EE_rdm
    '''
    
    
    
    


def Mutual(groundstate,i,j,L,dim):  #Quantum mutual information
        
     # creating "A" matrix
     I = sparse.csc_matrix(np.eye(dim))
     #%%
     ###################################
     #  Creating Ai  matrix
     ####################################
    
     A_i = dim*[[]]
     for d in range(dim):
         vec = np.zeros(dim)
         vec[d]=1
         if i!=0:
            S = vec
            for i_site in range(1,L):
                 if i_site==i:
                    s = I
                 else:
                    s = vec
                 S = sparse.kron(S,s, 'csc')
         else:
            S = I
            for i_site in range(L-1):
                S = sparse.kron(S,vec, 'csc')
         A_i[d] = S
     
     
     #%%
     ###################################
     #  Creating Aj  matrix
     ####################################
     A_j = dim*[[]]
     for d in range(dim):
         vec = np.zeros(dim)
         vec[d]=1
         if i!=0:
            S = vec
            for site in range(1,L):
                 if site==i:
                    s = I
                 else:
                    s = vec
                 S = sparse.kron(S,s, 'csc')
         else:
            S = I
            for site in range(L-1):
                S = sparse.kron(S,vec, 'csc')
         A_j[d] = S
        
     
     #%%
     ###################################
     #  Creating Aij  matrix
     ####################################
     A_ij = dim*[[]]
     for d in range(dim):
         
         vec = np.zeros(dim)
         vec[d]=1
         if i!=0 and j!=0:
            S = vec
            for site in range(1,L):
                if site==i or site==j:
                    s = I
                else:
                    s = vec
                S = sparse.kron(S,s, 'csc')
         else:
            S = I
            for site in range(1,L):
                if site==i or site==j:
                    s = I
                else:
                    s = vec
                S = sparse.kron(S,s, 'csc')
         A_ij[d] = S
        
         #print(A_ij[0].toarray())
         #print('===============')
         #for i in range(A_ij[0].shape[0]):
         #     aa = sparse.csc_matrix.getrow(A_ij[0],i)
         #     #print(aa)
         #     print(aa.nonzero()[1][0])
        
        
        
     # creating reduced density matrix
     rho_i = np.zeros((dim,dim),dtype=complex)
     for n,l in enumerate(A_i[0].nonzero()[1]):
         for m,k in enumerate(A_i[0].nonzero()[1]):
             rho_i[n,m] = groundstate[l].conj() * groundstate[k]
     for n,l in enumerate(A_i[1].nonzero()[1]):
         for m,k in enumerate(A_i[1].nonzero()[1]):
             rho_i[n,m] += groundstate[l].conj() * groundstate[k]
    
            
            
     rho_j = np.zeros((dim,dim),dtype=complex)
     for n,l in enumerate(A_j[0].nonzero()[1]):
         for m,k in enumerate(A_j[0].nonzero()[1]):
             rho_j[n,m] = groundstate[l].conj() * groundstate[k]
     for n,l in enumerate(A_j[1].nonzero()[1]):
         for m,k in enumerate(A_j[1].nonzero()[1]):
             rho_j[n,m] += groundstate[l].conj() * groundstate[k]
    
            
            
     rho_ij = np.zeros((2**dim,2**dim),dtype=complex)
     for n,l in enumerate(A_ij[0].nonzero()[1]):
          for m,k in enumerate(A_ij[0].nonzero()[1]):
             rho_ij[n,m] = groundstate[l].conj() * groundstate[k]
     for n,l in enumerate(A_ij[1].nonzero()[1]):
          for m,k in enumerate(A_ij[1].nonzero()[1]):
             rho_ij[n,m] += groundstate[l].conj() * groundstate[k]
            
            
     eig_i = sp.linalg.eigh(rho_i,eigvals_only=True)
     #Mu_i = -np.nansum(np.real(eig_i*np.log(eig_i)))
     eig_i[eig_i <= 0] = 1e-20
     svn  = -np.sum(np.real(eig_i*np.log(eig_i)))
     if math.isnan(svn):
         Mu_i = svn
     else:
         Mu_i = svn
    
    
     eig_j = sp.linalg.eigh(rho_j,eigvals_only=True)
     #Mu_j = -np.nansum(np.real(eig_j*np.log(eig_j)))
     eig_j[eig_j <= 0] = 1e-20
     svn  = -np.sum(np.real(eig_j*np.log(eig_j)))
     if math.isnan(svn):
         Mu_j = svn
     else:
         Mu_j = svn
         
     eig_ij = sp.linalg.eigh(rho_ij,eigvals_only=True)
     #Mu_ij = -np.nansum(np.real(eig_ij*np.log(eig_ij)))
     eig_ij[eig_ij <= 0] = 1e-20
     svn  = -np.sum(np.real(eig_ij*np.log(eig_ij)))
     if math.isnan(svn):
         Mu_ij = svn
     else:
         Mu_ij = svn
    
     return Mu_i + Mu_j - Mu_ij

    

#def Con(rho):
#
#    sy = np.array([[0, -1j], [1j, 0]])
#    rho_tilde = np.kron(sy,sy) @ rho @np.kron(sy,sy)
#    R = rho @ rho_tilde
#    eigs = np.linalg.eigvals(R)
#    eigs = np.sort(eigs.real)[::-1]
    #egv_max = max( eigs )
#
#    return max(eigs[0]-eigs[1]-eigs[2]-eigs[3],0)
