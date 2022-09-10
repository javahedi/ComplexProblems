lanczos_path = os.path.join(os.getcwd(),"lanczos")
sys.path.insert(0,lanczos_path)
#os.environ['KMP_DUPLICATE_LIB_OK']='True' # uncomment this line if omp error occurs on OSX for python 3
#os.environ['OMP_NUM_THREADS']='4' # set number of OpenMP threads to run in parallel
#os.environ['MKL_NUM_THREADS']='4' # set number of MKL threads to run in parallel


import scipy.sparse      as sparse
import matplotlib.pylab  as plt
import gen_diagprojector as GD
import gen_op_total      as GOT
import gen_operator      as GO
import gen_nn_int        as Gint
import numpy             as np
from multiprocessing     import Process
from multiprocessing     import current_process
from lanczos             import FTLM_static_iteration
from lanczos             import lanczos_full
from lanczos             import lanczos_iter
from lanczos             import lin_comb_Q_T
from lanczos             import expm_lanczos
from joblib              import Parallel
from joblib              import delayed
from numba               import jit
from numba               import njit
from numba               import prange
from scipy               import linalg
from time                import sleep
from time                import time
import sys
import os
import gc 






# Super (parent) class
class FockParafermion:
    '''
    *** Class for global parameters of the Fock Parafermions  ***
    '''
    # Class attribute
    class_name = "FockParafermions"
    def __init__(self,L=8,p=3,W=0.0,g=0.0,projection='True'):
        self.L          = L               #int(np.loadtxt('ChainLength.txt'))
        self.p          = p               #int(np.loadtxt('Flavour.txt'))
        self.projection = projection
        self.W          = W               # np.loadtxt('Wlist')   
        self.g          = g               # np.loadtxt('Glist')                                         
                            # generate  random on-site

    def operator(self):
        # %% Generate operators
        self.f0, self.fu, self.nop = GO.gen_operator(self.L,self.p)
        self.Ntot                  = GOT.gen_op_total(self.nop)
        if self.projection:
            val    = self.L//2    #np.loadtxt('ParticleNumber.txt')
            self.P = GD.gen_diagprojector(self.Ntot.diagonal(), val)

    def process_start(self):
        p = current_process()
        print('Starting:', p.name, p.pid)
        sys.stdout.flush()

    def process_finished(self):
        p = current_process()
        print('Finished:', p.name, p.pid)
        sys.stdout.flush()
    
    def Hamiltonian(self):
        self.operator()
        self.g_list     = [self.g] * self.L                                            
        self.m_list    = np.random.uniform(-self.W, self.W, self.L)
        self.Hint      = Gint.gen_nn_int(self.fu,self.nop,self.g_list,self.m_list,bc='obc')  # creat hamiltonian
        if self.projection:
            self.Hint  = self.P @ self.Hint @ self.P.T                             # projecting Hmailtonian
    
    def diagonaization(self,full='True'):
        self.Hamiltonian()
        if full:
            self.energy, self.evecs = linalg.eigh(self.Hint.todense())             # full diagonaization of H


    def __str__(self):
        return f"\
                \nPareametrs are set as follows:    \
                \n==============================    \
                \nLenght        : {self.L}          \
                \nDisorder      : {self.W}          \
                \nG             : {self.g}          \
                \nProjection    : {self.projection} \
                \nHilber dim    : {self.f0[0].shape}\
                \nProjected dim : {self.P.shape}    \
                \nElapsed time  : {self.elapsed}  seconds"
    

# Child class
class transport(FockParafermion):

    # Instance method
    def autocorr_ED(self,tmax=101):
        self.process_start()
        self.diagonaization()
        timeList = np.linspace(0.,tmax,201,endpoint=True)
        start    = time()
        # particle operator
        Nl       = self.nop[self.L//2]
        Nl       = self.P @ Nl @ self.P.T
        cc       = self.evecs.conj().T @ Nl @ self.evecs
        try :
            # Cython fast loop
            import fast_loop
            autoco = fast_loop.TimeDomian(timeList,self.energy,cc)
        except ModuleNotFoundError:
            print(" program running on NUMBA ")
            autoco = numpa_time(timeList,self.energy,cc,cc)
        autoco  /=self.energy.shape[0]
        np.savetxt('autocorr_ed.txt',autoco.real)
        self.elapsed =  time()-start
        self.process_finished()



    # Instance method
    def autocorr_QDT(self,tmax=101,method='rk4'):
        self.process_start()
        self.diagonaization()
        start    = time()
        timeList = np.linspace(0.,tmax,201,endpoint=True)
        dt       = timeList[1] - timeList[0]
        Nl =  self.nop[self.L//2]
        Nl =  self.P @ Nl @ self.P.T

        mu,sigma=0,1
        a    = np.random.normal(mu,sigma,size=self.f0[0].shape[0])
        b    = np.random.normal(mu,sigma,size=self.f0[0].shape[0])
        phi  = (a  + 1j*b)
        phi  = self.P @ phi
        phi /= np.linalg.norm(phi)
        #phi  = np.sqrt(Nl) @ phi/np.sqrt(np.linalg.norm(phi))
        autoco = []
        for _ in timeList:
            if method=='rk4':
                phi = RK4(phi,self.Hint,dt)
                autoco.append((phi.conj().T @ Nl @ phi))
            if method=='Krylov':
                E, V, Q_T = lanczos_iter(self.Hint,phi,20 )
                phi       = expm_lanczos(E, V, Q_T, a=-1j*dt)
                autoco.append( (phi.conj().T @ Nl @ phi )/2 )
        np.savetxt('autocorr_qdt.txt',np.array(autoco).real)   
        self.elapsed =  time()-start
        self.process_finished()


    # Instance method
    def mean_square_distance(self,tmax=101):
        self.process_start()
        self.diagonaization()
        timeList = np.logspace(0,1.7,tmax,base=10)
        start    = time()
        msd = 0
        for s in range(self.L):
            Ns =  self.P @ self.nop[s] @ self.P.T
            cs =  self.evecs.conj().T @ Ns @ self.evecs
            for r in range(self.L):
                Nr =  self.P @ self.nop[r] @ self.P.T
                cr =  self.evecs.conj().T @ Nr @ self.evecs
                try :
                    #Cython fast loop
                    import fast_loop
                    ctr_timeT = fast_loop.TimeDomian2(timeList,self.energy,cr,cs)
                except ModuleNotFoundError:
                    print(" program running on NUMBA ")
                    ctr_timeT = numpa_time(timeList,self.energy,cr,cs) 
                msd += (s-r)**2 * ctr_timeT
        msd /=  self.energy.shape[0]
        np.savetxt('mean_square_distance_sample.txt',msd.real/self.L)
        self.elapsed =  time()-start
        self.process_finished()




    # Instance method
    def weighted_correlation(self,tmax=101):
        
        self.process_start()
        self.diagonaization()
        start    = time()
        timeList = np.logspace(0,1.7,tmax,base=10)
        N0 =  self.P @ self.nop[0] @ self.P.T
        c0 =  self.evecs.conj().T @ N0 @ self.evecs
        x2 = 0
        for r in range(1,self.L):
            Nr =  self.P @ self.nop[r] @ self.P.T
            cr =  self.evecs.conj().T @ Nr @ self.evecs
            Cr0_time0 = np.sum(self.evecs.conj().T @ Nr @ N0 @ self.evecs) # at time zero
            try :
                #Cython fast loop
                import fast_loop
                Cr0_timeT = fast_loop.TimeDomian2(timeList,self.energy,cr,c0)
            except ModuleNotFoundError:
                print(" program running on NUMBA ")
                Cr0_timeT = numpa_time(timeList,self.energy,cr,c0)
            x2 += r**2 * (Cr0_timeT-Cr0_time0*np.ones_like(timeList))
        x2 /=self.energy.shape[0]
        self.elapsed =  time()-start
        np.savetxt('weighted_correlation_sample.txt',x2.real)
        self.process_finished()


    # Instance method
    def conductivity_ED(self,omega_max=100):
        self.process_start()
        self.diagonaization()
        start  = time()
        Omega  = np.logspace(-2,0.5,omega_max,base=10)
        # for OBC = L-1,  PBC = L
        bound   = [(site, (site+1)%self.L) for site in range(self.L-1)]
        current = sparse.csr_matrix(self.fu[0].shape,dtype=np.complex64)
        for ind,val in enumerate(bound):
            current += -1j * (1-self.g_list[ind]) * (self.fu[val[0]]   @ self.fu[val[1]].T )
            current +=  1j * (1-self.g_list[ind]) * (self.fu[val[0]].T @ self.fu[val[1]]   )
            current += -1j *  2*self.g_list[ind]  * (self.fu[val[0]]   @ self.fu[val[0]]   @  self.fu[val[1]].T @  self.fu[val[1]].T )
            current +=  1j *  2*self.g_list[ind]  * (self.fu[val[0]].T @ self.fu[val[0]].T @  self.fu[val[1]]   @  self.fu[val[1]]   )

        current = self.P              @ current @ self.P.T
        cc      = self.evecs.conj().T @ current @ self.evecs

        # sigma is brodening 
        sigma  = self.W * np.sqrt(self.L)/self.energy.shape[0]/100
        sigma2 = np.diff(Omega)
        sigma2 = np.hstack((sigma2,sigma2[-1]))
        try:
            import fast_loop
            con_ed,con_ed2 = fast_loop.OmegaDomian2(Omega,self.energy,cc,sigma,sigma2)
        except ModuleNotFoundError:
            print(" program running on NUMBA ")
            con_ed,con_ed2 = numpa_omega(Omega,self.energy,cc,sigma,sigma2)
        con_ed  = con_ed /self.energy.shape[0]/self.L
        con_ed2 = con_ed2/self.energy.shape[0]/self.L
        np.savetxt('conductivity_ed.txt', con_ed.real)
        np.savetxt('conductivity2_ed.txt',con_ed2.real)
        self.elapsed =  time()-start
        self.process_finished()


    # Instance method
    def conductivity_QDT(self,omega_max=100,method='rk4'):
        self.process_start()
        self.diagonaization()
        start  = time()
        
        Omega  = np.logspace(-2,0.5,omega_max,base=10)
        # for OBC = L-1,  PBC = L
        bound   = [(site, (site+1)%self.L) for site in range(self.L-1)]
        current = sparse.csr_matrix(self.fu[0].shape,dtype=np.complex64)
        for ind,val in enumerate(bound):
            current += -1j * (1-self.g_list[ind]) * (self.fu[val[0]]   @ self.fu[val[1]].T )
            current +=  1j * (1-self.g_list[ind]) * (self.fu[val[0]].T @ self.fu[val[1]]   )
            current += -1j *  2*self.g_list[ind]  * (self.fu[val[0]]   @ self.fu[val[0]]   @  self.fu[val[1]].T @  self.fu[val[1]].T )
            current +=  1j *  2*self.g_list[ind]  * (self.fu[val[0]].T @ self.fu[val[0]].T @  self.fu[val[1]]   @  self.fu[val[1]]   )

        current = self.P @ current @ self.P.T


        mu,sigma=0,1
        a    = np.random.normal(mu,sigma,size=self.f0[0].shape[0])
        b    = np.random.normal(mu,sigma,size=self.f0[0].shape[0])
        phi  = (a  + 1j*b)
        phi  = self.P @ phi
        phi /= np.linalg.norm(phi)
        left  = phi
        right = current @ phi

        tmax = 100
        timeList = np.linspace(0.,tmax,2000,endpoint=False)
        dt = timeList[1] - timeList[0]
        curent_curent = []
        for _ in timeList:
            if method=='rk4':
                left  = RK4(left, self.Hint,-1j*dt)
                right = RK4(right,self.Hint,-1j*dt)
                curent_curent.append((left.conj().T @ current @ right)/self.L)
            if method=='Krylov':
                # left
                E, V, Q_T = lanczos_iter(self.Hint,left,20)
                left = expm_lanczos(E, V, Q_T, a = -1j*dt )
                # right
                E, V, Q_T = lanczos_iter(self.Hint,right,20)
                right = expm_lanczos(E, V, Q_T, a = -1j*dt )
                curent_curent.append( (left.conj().T @ current @ right)/self.L) 

        curent_curent = np.array(curent_curent)
        con_dqt  = np.zeros(Omega.shape,dtype=np.float64)
        for id, w in enumerate(Omega):
            Current  = np.cumsum(curent_curent* np.exp(1j*w*timeList) * dt)[-1]
            con_dqt[id] = Current.real

        np.savetxt('conductivity_qdt.txt',np.array(con_dqt).real)   
        self.elapsed =  time()-start
        self.process_finished()


##########################
##########################
##########################
@njit(parallel=True, fastmath=True)
def numpa_time(timeList,ee,vec1,vec2):
    out = np.zeros(len(timeList),dtype=np.complex64)
    for i in prange(len(ee)):
        for j in prange(len(ee)):
            out +=  np.exp(1j*(ee[j]-ee[i])*timeList) * vec1[i,j] * vec2[j,i]
    return out   
##########################
##########################
##########################
@njit(parallel=True, fastmath=True)
def numpa_omega(Omega,ee,vec,s1,s2):
    out  = np.zeros(len(Omega),dtype=np.float64)
    out2 = np.zeros(len(Omega),dtype=np.float64)
    for i in prange(len(ee)):
        for j in prange(len(ee)):
            out  += LorentzSmearing(Omega, ee[j]-ee[i],s1) * abs(vec[i,j])**2
            out2 += LorentzSmearing(Omega, ee[j]-ee[i],s2) * abs(vec[i,j])**2
    return out, out2  
##########################
##########################
##########################
@njit
def LorentzSmearing(x,x0,sigma):
    return sigma /((x - x0)**2 + sigma**2)/np.pi
##########################
##########################
##########################
@njit
def RK4(fk,ham,dt):
    out = fk
    for k in range(4):
        fk *= dt*ham/(k+1)
        out +=fk
    return out
################################################
################################################

if __name__ == "__main__":
    print(FockParafermion.__doc__)
    model = transport(6,3,1.0,0.0,'True')
    
    jobs = []
    
    ##################
    job1 = Process(name="autocorr_ED", target=model.autocorr_ED())
    jobs.append(job1)
    job1.daemon = True
    job1.start()
    
    ##################
    #job2 = Process(name="mean_square_distance", target=model.mean_square_distance())
    #jobs.append(job2)
    #job2.daemon = True
    #job2.start()
    
    ##################
    job3 = Process(name="conductivity_ED",      target=model.conductivity_ED())
    jobs.append(job3)
    job3.daemon = True
    job3.start()
    
    
    job1.join()
    #job2.join()
    job3.join()

    

    

    print ('job1 is alive :', job1.is_alive())
    #print ('job2 is alive :', job2.is_alive())
    print ('job3 is alive :', job3.is_alive())


    
    



    #print(model)
