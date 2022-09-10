
import numpy as np
from scipy.sparse import csc_matrix
import matplotlib.pyplot as plt
from numba import njit

omega = np.exp(1j*2*np.pi/3)

#@njit
def create_operator(a,localdim=None,n=None,name=None):
    """Create the different operators"""
    operators = dict() # empty dictionary
    if name     is None: raise
    if n        is None: raise
    if localdim is None: raise
    ids = [np.identity(j,dtype=np.complex) for j in localdim] # identities

    # create operators in  for all sites
    for i in range(n):
        if a.shape[0]!=localdim[i]: raise
        op = one2many(ids,a,i)                           # one to many body
        operators[(name,i)] = op                    # store in the dictionary
    return operators

 
#@njit
def one2many(ids,op=None,i=-1):
    """Function to transform to many body basis given identity operators"""
    tmp = np.zeros((1,1),dtype=np.complex) # initialize
    tmp[0,0] = 1.0
    for j in range(len(ids)): # loop over sites
        if i!=j:
            op2 = ids[j] # identity
        else:
            op2 =  op    # operator
        tmp = np.kron(tmp,op2) # tensor product
    tmp = csc_matrix(tmp) # return operator
    tmp.eliminate_zeros()
    return tmp

#@njit
def Parafermion_Chain(n,Z):
    # n : number of sites
    # create all operators
    operators = dict() # empty dictionary
    localdim = [Z for i in range(n)]
    zero = np.zeros((Z,Z),dtype=np.complex) # get zero

    #%%
    N = zero.copy()
    for i in range(Z): N[i,i] = i
    #%%
    S = zero.copy()
    for i in range(Z-1):
        S[i,i+1] = 1.0
    S[Z-1,0] = 1.0
    #%%
    T = zero.copy()
    for i in range(Z): T[i,i] = np.exp(i*1j*2*np.pi/Z)
    Td = np.conjugate(T.T)

    Ntot   = create_operator(N  ,localdim, n, name="N"     )
    Sig    = create_operator(S  ,localdim, n, name="Sig"   )
    SigDag = create_operator(S.T,localdim, n, name="SigDag")
    Tau    = create_operator(T  ,localdim, n, name="Tau"   )
    TauDag = create_operator(Td ,localdim, n, name="TauDag")

    return Ntot, Sig, SigDag, Tau, TauDag

#@njit
def test_commutation(OP1,OP2=[],name=[],Z=3):
    """Perform a test of the commutation relations"""
    if OP2==[]: OP2=OP1
    if name==[]: raise
        
#    Id = self.get_operator(self.Id)
    n = len(OP1) # number of sites
    ntries = 8 # number of tries
    omega = np.exp(1j*2*np.pi/Z)
    for _ in range(ntries):
        i = np.random.randint(n)
        j = np.random.randint(n)
        if i>=j: continue
        d = OP1[(name[0], i)]@OP2[(name[1], j)] - omega*OP2[(name[1], j)]@OP1[(name[0], i)]
        
        if np.max(np.abs(d))>1e-7:
            print(f"{name[0]}, {name[1]} failed at sites : ({i},{j})")
            print(f"===> {np.max(np.abs(d))}")
            raise
    print("Commutation test passed")



if __name__ == "__main__":
    n,Z=3,3
    Ntot, Sig, SigDag, Tau, TauDag = Parafermion_Chain(n,Z)
    
    #%% plot
    fig, (ax0, ax1) = plt.subplots(figsize=(10, 3), ncols=2)
    fu0 = ax0.imshow(Sig[('Sig', 0)].todense().real, cmap='Blues', interpolation='none')
    fu1 = ax1.imshow(Sig[('Sig', 1)].todense().real, cmap='Blues', interpolation='none')
    fig.colorbar(fu0, ax=ax0)
    fig.colorbar(fu1, ax=ax1)
    #ax0.spy(Sig[('Sig', 0)].todense().real)
    #ax1.spy(SigDag[('SigDag', 0)].todense().real)
    plt.show()
    
    #%% test
    test_commutation(Sig,Sig,name=['Sig','Sig'])
    test_commutation(Tau,Sig,name=['Tau','Tau'])
    test_commutation(Sig,Tau,name=['Sig','Tau'])
