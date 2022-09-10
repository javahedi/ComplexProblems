import numpy as np
#import math
from scipy.linalg import eigh
import Hamiltonian
import itertools
import os

def SelfCon(Dim,U,Density_up,Density_dn,Nx,Ny,beta):
    itr=0
    nocc =  Dim
    #filling = int(round(Nx*Ny*(nocc+2)))    ##  6/8=2*(3/4)
    filling = int(Nx*Ny*nocc+0)   ## half filling  4/8=2*(1/4)

    changeUp=1
    changeDn=1
    #Reciproque lattice: in Unit of 2pi/a  a=2.456 Angstrom
    b1 = np.array([6.91756085e-02,  -3.68663557e-02])
    b2 = np.array([-2.66060024e-03,   7.83410147e-02])
   
    while changeUp > 1e-6 and changeDn > 1e-6:
    
        itr += 1 
        New_Density_up = np.zeros(Dim) 
        New_Density_dn = np.zeros(Dim) 
        
        #################################################
        #################################################
        #################################################
        EnUp = []
        EnDn = []
        s = 0
        for nx, ny in itertools.product(range(Nx), range(Ny)):
               
            K = (nx/Nx)*b1+(ny/Nx)*b2
            HK_up, HK_dn = Hamiltonian.Hamil(Dim,U,Density_up,Density_dn,K,s,0)
            
            en = eigh(HK_up,  eigvals_only=True)            
            EnUp.extend(en)
            en = eigh(HK_dn,  eigvals_only=True)
            EnDn.extend(en)
            
            s +=1
            
        eigud = np.sort(EnUp+EnDn)
        fermi = (eigud[filling] + eigud[filling-1])/2
        ####################################################
        ####################################################
        ####################################################
        s = 0
        gap_up = []
        gap_dn = []
        E_GState = 0
        for nx, ny in itertools.product(range(Nx), range(Ny)):

            K = (nx/Nx)*b1+(ny/Nx)*b2            
            HK_up, HK_dn= Hamiltonian.Hamil(Dim,U,Density_up,Density_dn,K,s,1)
            
            eigs_u, vecs_u = np.linalg.eigh(HK_up)
            eigs_d, vecs_d = np.linalg.eigh(HK_dn)

            
            Su,Sd = 0,0
            for j in range(Dim):
                if eigs_u[j]< fermi:
                    E_GState += eigs_u[j]
                    Su +=1
                    #for i in range(Dim):
                    New_Density_up[:] +=  abs(vecs_u[:,j])**2# * vecs_u[:,j].conj()
                if eigs_d[j]< fermi:
                    E_GState += eigs_d[j]
                    Sd +=1
                    #for i in range(Dim):
                    New_Density_dn[:] +=  abs(vecs_d[:,j])**2# * vecs_d[i,j].conj()
        
                                   
            gap_up.append(eigs_u[Su]-eigs_u[Su-1])
            gap_dn.append(eigs_d[Sd]-eigs_d[Sd-1])
            s +=1       

        
        ################################################
        #######      End of FOR loop     #############
        ################################################
   
        ################################################
        ####### Finding New Density          ###########
        ################################################
        p = 0.7 # mixing parameter
        New_Density_up = p*Density_up + (1-p)*New_Density_up/(Nx*Ny)
        New_Density_dn = p*Density_dn + (1-p)*New_Density_dn/(Nx*Ny)
          
       
        changeUp = max(abs(New_Density_up-Density_up))  
        changeDn = max(abs(New_Density_dn-Density_dn)) 
        
        Density_up=New_Density_up
        Density_dn=New_Density_dn

       
        print ('Ite = ' + str(itr) + ':' +' chU = ' + str(changeUp) +\
               ', chD = ' + str(changeDn) + '',flush=True)
        
        os.system("rm HK_*")
    
    E_GState = E_GState/(Nx*Ny)
    #E_GState -= U*(New_Density_up*New_Density_dn).sum()
    return itr, Density_up, Density_dn, E_GState/Dim, min(gap_up), min(gap_dn)
        
