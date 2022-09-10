import numpy as np
import Hamiltonian
#from scipy.linalg import eigh
#import GenerateKpoint
import Structure
import itertools
import time


def SelfCon(Dim,U,Density_up,Density_dn,Nx,Ny):
    itr=0
    
    changeUp=1
    changeDn=1
    #dEG = 1
    ###############################################
    ##############################################
    m,n=8,9
    t0=-0.75 
    
    filename,L1,L2,G2,G1 = Structure.UnitCell(m,n)

    H00,HL1,HL2,HL1L2,HL1_L2 =\
        Hamiltonian.Hamil(t0,Dim,filename,m,n,L1,L2)
    ###############################################
    ###############################################
    Nk = Nx*Ny
    Eu = np.zeros((Nk,Dim))
    Ed = np.zeros((Nk,Dim))
    Vu = np.zeros((Nk,Dim,Dim),complex)
    Vd = np.zeros((Nk,Dim,Dim),complex)
    H00_u , H00_d = H00, H00
    
    
    
    '''
    np.random.seed(312345)
    G = np.array([0.0,  0.0])
    M =  (G1+G2)/2
    K1 = ((1/3)*G1+(2/3)*G2)
    K2 = [-K1[0],K1[-1]]
    K3 = -((2/3)*G1+(1/3)*G2)
    K4 = [K2[0],-K2[-1]]
    K5 = [K1[0],-K1[-1]]
    K6 = - K3

    kx = [G[0]]
    ky = [G[1]]
    
    for i in range(nk):
        
        x,y=GenerateKpoint.trisample(G, K1, K2)
        kx.append(x)
        ky.append(y)
    
        #x,y=GenerateKpoint.trisample(G, K2, K3)
        #kx.append(x)
       	#ky.append(y)
           
        #x,y=GenerateKpoint.trisample(G, K3, K4)
        #kx.append(x)
        #ky.append(y)
        
        #x,y=GenerateKpoint.trisample(G, K4, K5)
        #kx.append(x)
        #ky.append(y)
        
        #x,y=GenerateKpoint.trisample(G, K5, K6)
        #kx.append(x)
        #ky.append(y)
        
        #x,y=GenerateKpoint.trisample(G, K6, K1)
        #kx.append(x)
        #ky.append(y)
    ###########################################
    ###########################################
    '''    
    while changeUp > 1e-6 and changeDn > 1e-6:
    
        start = time.time()
        itr += 1 
        New_Density_up = np.zeros(Dim) 
        New_Density_dn = np.zeros(Dim)

       
        Etot=[]
        k=0
        for nx, ny in itertools.product(range(Nx), range(Ny)):
                      
            K = (nx/Nx)*G1+(ny/Ny)*G2
            KL1,  KL2    = np.dot(K,L1),    np.dot(K,L2)
            KL1L2,KL1_L2 = np.dot(K,L1+L2), np.dot(K,L1-L2)

            ################################################
            #######            H(k)              ###########
            ################################################
            np.fill_diagonal(H00_u, U*Density_dn)
             
            H_up = H00_u+HL1.T*np.exp(-1j*KL1)+HL1*np.exp(1j*KL1)+\
              	      HL2.T*np.exp(-1j*KL2)+HL2*np.exp(1j*KL2)+\
                      HL1L2.T*np.exp(-1j*KL1L2)+HL1L2*np.exp(1j*KL1L2)+\
                      HL1_L2.T*np.exp(-1j*KL1_L2)+HL1_L2*np.exp(1j*KL1_L2)
                      
            en_u, vec_u=np.linalg.eigh(H_up)
            idx = en_u.argsort()[::]   
            Eu[k,:] = en_u[idx]
            Vu[k,:,:] = vec_u[:,idx] 
            
            np.fill_diagonal(H00_d, U*Density_up)             
            H_dn = H00_d+HL1.T*np.exp(-1j*KL1)+HL1*np.exp(1j*KL1)+\
              	      HL2.T*np.exp(-1j*KL2)+HL2*np.exp(1j*KL2)+\
                      HL1L2.T*np.exp(-1j*KL1L2)+HL1L2*np.exp(1j*KL1L2)+\
                      HL1_L2.T*np.exp(-1j*KL1_L2)+HL1_L2*np.exp(1j*KL1_L2)
            
            en_d, vec_d=np.linalg.eigh(H_dn)
            idx = en_d.argsort()[::]              
            Ed[k,:] = en_d[idx]
            Vd[k,:,:] = vec_d[:,idx]
            
            Etot.extend(en_u)
            Etot.extend(en_d)
            
            k +=1
              

        Etot.sort()
        fermi=(Etot[round(Dim*Nk)]+Etot[round(Dim*Nk)-1])/2            
           
        ####################################################
        ####################################################
        ####################################################
        gap_up = []
        gap_dn = []
        E_GState = 0        
        for k in range(Nk):                   
            Su,Sd = 0,0
            for j in range(Dim):
                if Eu[k,j]< fermi:
                    E_GState += Eu[k,j]
                    Su +=1   
                    for i in range(Dim):
                        New_Density_up[i] +=  abs(Vu[k,i,j])**2
                if Ed[k,j]< fermi:
                    E_GState += Ed[k,j]
                    Sd +=1                    
                    for i in range(Dim):
                        New_Density_dn[i] +=  abs(Vd[k,i,j])**2
    
            gap_up.append(Eu[k,Su]-Eu[k,Su-1])
            gap_dn.append(Ed[k,Sd]-Ed[k,Sd-1])               
                           
        ################################################
        #######      End of FOR Loop     #############
        ################################################
   
        ################################################
        ####### Finding New Density          ###########
        ################################################
        
        
        p=0.7
        New_Density_up=p*Density_up+(1-p)*New_Density_up/Nk
        New_Density_dn=p*Density_dn+(1-p)*New_Density_dn/Nk         
        
        changeUp=max(abs(New_Density_up-Density_up))
        changeDn=max(abs(New_Density_dn-Density_dn)) 
        
        Density_up = New_Density_up 
        Density_dn = New_Density_dn 
        
        Charge = Density_up+Density_dn
        criteria = 1-sum(Charge)/Dim
        print('itr = '+str(itr)+'  criteria = '+str(criteria)+\
              '  changeUp = '+str(changeUp)+'  changeDn = '+str(changeDn)+''\
              '  Time = '+str(time.time()-start)+'')
        
        
    ################################################
    #######      End of while loop     #############
    ################################################
    E_GState = E_GState/Nk/Sd/Su
    E_GState -= U*(Density_up*Density_up).sum()/Dim     

    return itr, Density_up, Density_dn, E_GState, min(gap_up), min(gap_dn)

        
