import numpy as np
#import math
import Hamiltonian
#import itertools

def SelfCon(t,U,Density_up,Density_down,kxlist,kylist,beta):
    itr=0;
    Dim = len(Density_down)
    changeUp=1;
    changeDown=1;
    while ((changeUp > 1e-6) & (changeDown > 1e-6)):
        #print(itr, changeUp)
        itr += 1 
        New_Density_up = np.zeros(Dim) 
        New_Density_down = np.zeros(Dim)
 
        for kx,ky in zip(range(len(kxlist)), range(len(kylist))):
        #for kx, ky in itertools.product(range(len(kxlist)), repeat = 2):
                k = [kx,ky]
                eigs_u, vecs_u, eigs_d, vecs_d = \
                Hamiltonian.Hamil(t,U,Density_up,Density_down,k)
                for j in range(Dim):
                  prefac_u = 0.5*(1-np.tanh(beta*(eigs_u[j])/2))
                  prefac_d = 0.5*(1-np.tanh(beta*(eigs_d[j])/2))
                  for i in range(Dim):
                        New_Density_up[i] +=  prefac_u *vecs_u[i][j]*vecs_u[i][j].conj()
                        New_Density_down[i] += prefac_d *vecs_d[i][j]*vecs_d[i][j].conj()
           
                
        ################################################
        #######      End of FOR loop     #############
        ################################################
   
        ################################################
        ####### Finding New Density          ###########
        ################################################
        p=0.7 # mixing parameter
        New_Density_up=p*Density_up+(1-p)*New_Density_up/len(kxlist)
        New_Density_down=p*Density_down+(1-p)*New_Density_down/len(kxlist)
          
        changeUp=np.sum(np.absolute(np.subtract(New_Density_up,Density_up)))/Dim    
        changeDown=np.sum(np.absolute(np.subtract(New_Density_down,Density_down)))/Dim
          
        Density_up=New_Density_up
        Density_down=New_Density_down
        
        


    return itr,Density_up, Density_down
        
