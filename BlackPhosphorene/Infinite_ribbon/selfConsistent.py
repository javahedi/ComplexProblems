import numpy as np

#from scipy.linalg import eigh


def SelfCon(Dim,U,Dens_up,Dens_dn,H00,H10,Nk,pp):
    
    H00_u,H00_d=H00,H00
    Eu=np.zeros((Nk,Dim))
    Ed=np.zeros((Nk,Dim))
    Vu=np.zeros((Nk,Dim,Dim),complex)
    Vd=np.zeros((Nk,Dim,Dim),complex)
    b = 3.32 * pp
    ky=np.linspace(-np.pi/b,np.pi/b,Nk)
    beta = 1e5
    itr=0
    changeUp=1
    changeDn=1
    while changeUp > 1e-6 and changeDn > 1e-6:
        itr = itr +1
        New_Dens_up = np.zeros(Dim) 
        New_Dens_dn = np.zeros(Dim)
        ################################################
        #######   First diagonalizing over the FBZ   ###
        #######   to find Femrmi in whole FBZ        ### 
        ################################################
        Etot=[]
        for k in range(Nk):  
            
            np.fill_diagonal(H00_u, U*(Dens_dn-1/2-Dens_dn*Dens_up))
             
            HH_u = H00_u+np.transpose(H10).conj()*np.exp(-1j*ky[k]*b)+\
                   H10*np.exp(1j*ky[k]*b)
            en_u, vec_u=np.linalg.eigh(HH_u)
            idx = en_u.argsort()[::]
            Eu[k,:] = en_u[idx]
            Vu[k,:,:] = vec_u[:,idx] 
            
            
            np.fill_diagonal(H00_d, U*(Dens_up-1/2-Dens_dn*Dens_up))
            HH_d = H00_d+np.transpose(H10).conj()*np.exp(-1j*ky[k]*b)+\
                   H10*np.exp(1j*ky[k]*b)
            en_d, vec_d=np.linalg.eigh(HH_d)
            idx = en_d.argsort()[::]
            Ed[k,:] = en_d[idx]
            Vd[k,:,:] = vec_d[:,idx]
                
            Etot.extend(en_u)
            Etot.extend(en_d)
            
        Etot.sort()
        fermi=(Etot[round(Dim*Nk)]+Etot[round(Dim*Nk)-1])/2

        ################################################
        ####### Finding New Dens          ###########
        ################################################
        
        E_GState = 0
        New_dens_up = np.zeros(Dim)
        New_dens_dn = np.zeros(Dim)
        
        '''
        for k in range(Nk):
            for j in range(Dim):
              prefac_u = 0.5 * (1-np.tanh(beta*(Eu[k,j]-fermi)/2))
              prefac_d = 0.5 * (1-np.tanh(beta*(Ed[k,j]-fermi)/2))
              E_GState += Eu[k,j] * prefac_u
              E_GState += Ed[k,j] * prefac_d
              New_dens_up[:] += prefac_u * abs(Vu[k,:,j])**2#*Vu[k,i,j].conj()
              New_dens_dn[:] += prefac_d * abs(Vd[k,:,j])**2#*Vd[k,i,j].conj()
        
        '''
 
        
        for k in range(Nk):
            for j in range(Dim):
                if Eu[k,j]< fermi:
                    E_GState += Eu[k,j]
                    New_Dens_up[:] +=  abs(Vu[k,:,j])**2
                if Ed[k,j]< fermi:
                    E_GState += Ed[k,j]
                    New_Dens_dn[:] +=  abs(Vd[k,:,j])**2
        
        ################################################
        #######      End of FOR Loop     #############
        ################################################
   
        ################################################
        ####### Finding New Dens          ###########
        ################################################
        p=0.7
        New_Dens_up=p*Dens_up+(1-p)*New_Dens_up/Nk
        New_Dens_dn=p*Dens_dn+(1-p)*New_Dens_dn/Nk
        
        changeUp=max(abs(New_Dens_up-Dens_up))
        changeDn=max(abs(New_Dens_dn-Dens_dn)) 
        
        Dens_up = New_Dens_up 
        Dens_dn = New_Dens_dn 
        
        #Charge = Dens_up+Dens_dn
        #criteria = 1-sum(Charge)/Dim
        print('U = '+str(U)+'  itr = '+str(itr)+'  changeUp = '+str(changeUp)+'')

    
    ################################################
    #######      End of while loop     #############
    ################################################
    

    return Dens_up, Dens_dn, E_GState/Nk/Dim
