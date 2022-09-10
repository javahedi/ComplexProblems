#from math import sqrt, pi, exp, cos, sin,log
import numpy as np
import os

##############################
def Hamil(Dim,U,Density_up,Density_dn,K,s,ss):
    
    HK_up = np.zeros((Dim,Dim),complex)
    HK_dn = np.zeros((Dim,Dim),complex)

    vec_k = [K[0] , K[1] , 0.0] # K
    
    f1 = open("KVECTOR", "w")
    f1.write(' '.join(str(item) for item in vec_k))
    f1.close
    
    f1 = open("KVECTOR", "r") 
    a = f1.read()
    f1.close
    


    # writing of Hk
    if ss==0:
       os.system('./writeHk.run > writeHk.out')
       os.system('mv HK HK_'+str(s)+'')
    site_i, site_j, real_HK, imag_HK = np.loadtxt('HK_'+str(s)+'', unpack=True)
    
    os.system("rm KVECTOR")
    ################################################
    #######  Call Up and Down  Hamiltonian  ########
    ################################################
    for i, j, c, d in zip(site_i, site_j, real_HK, imag_HK):
        HK_up[int(i)-1,int(j)-1] = c+1j*d
        HK_up[int(j)-1,int(i)-1] = c-1j*d
        
        HK_dn[int(i)-1,int(j)-1] = c+1j*d
        HK_dn[int(j)-1,int(i)-1] = c-1j*d
    
   
    np.fill_diagonal(HK_up, U*Density_dn)
    np.fill_diagonal(HK_dn, U*Density_up)
      
    return HK_up, HK_dn
