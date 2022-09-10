"""
Authors: Thu Phung + Javad Vahedi
Date : December 13, 2017  
Project:  The Hubbard model on Graphene sheet in k_space
"""
import numpy as np
from math import pi
import matplotlib.pyplot as plt
#from scipy import optimize
import time
import SelfConsistent
import Strain


##########################################
##########################################
##########################################
##########################################
#       Set parameters here

# Hopping parameters
t1 = -1.22  # eV
t2 = 3.665
t3 = -0.205
t4 = -0.105
t5 = -0.055
delta = 0.0  # gate
e0 = delta*abs(t1)
beta = 1e5   # temperature

# strain parameter
ex, ey, ez = 0.0, 0.0, 0.0
Hopping = [t1, t2, t3, t4, t5]
Hopping = Strain.strain(ex, ey, ez, Hopping)

mesh = 1000
ax = 4.45       
ay = 3.32
kxlist = np.linspace(-pi/ax,pi/ax,mesh)
kylist = np.linspace(-pi/ay,pi/ay,mesh)

 
Magnetization_a = []
Magnetization_b = []
E_GState = []
E_Gap_u = []
E_Gap_d = []

Urange = np.arange(2,3,0.1)
f1 = open('Magnetization_a.dat', 'w')
f2 = open('E_Gap_u.dat', 'w')
#f3 = open('E_GState.dat', 'w')


for U in Urange:
    ################################################
    #######         initial value           ########
    ################################################ 
    Density_up   = np.array([0.5, 0.0])
    Density_down = np.array([0.0, 0.5])
    
    start = time.clock()    
    ################################################
    #######  Call self-consistent           ########
    ################################################
    itr, Density_up, Density_down = \
                SelfConsistent.SelfCon(Hopping,U,Density_up,Density_down,kxlist,kylist,beta)

    
    ################################################
    #######      Single-Site Gap         ###########
    ################################################    
    E_Gap_u.append(np.absolute(U*Density_down[0]-U*Density_down[1]))    
    E_Gap_d.append(np.absolute(U*Density_up[0]-U*Density_up[1]))

    ################################################
    #######      Ground state Energy     ###########
    ################################################
#    Eu = 0
#    Ed = 0
#    for kx,ky in itertools.product(range(mesh), repeat = 2):
#    #for kx,ky in zip(kxlist,kylist):
#        k = [kxlist[kx], kylist[ky]]
#        eigs_u, vecs_u, eigs_d, vecs_d = \
#                Hamiltonian.Hamil(Hopping,U,Density_up,Density_down,k)
#        Eu+=eigs_u[0]
#        Ed+=eigs_d[0]
#        
#       
#    E_GState.append(((Eu+Ed)/mesh**2)-U*np.dot(Density_up,Density_down))

    
    ################################################
    #######      On-Site Magnetization   ###########
    ################################################  
    Sz_a=np.subtract(Density_up[0],Density_down[0])/2
    Magnetization_a.append(Sz_a)
  

    Sz_b=np.subtract(Density_up[1],Density_down[1])/2; 
    Magnetization_b.append(Sz_b)
    
   

    f1.write("%f  %f\n"  %(U, Magnetization_a[-1] ))
    f2.write("%f  %f\n"  %(U, E_Gap_u[-1] ))   
    #f3.write("%f  %f\n"  %(U, E_GState[-1] )) 
    
   
    print ('U = '+ str(U) +   ' ,itr = ' + str(itr) +\
           ', Magnetization_a = ' + str(Magnetization_a[-1]) +', Time = ' + str(time.clock() - start))
   
f1.close()
f2.close()
#f3.close()



plt.subplot(1,1,1)
plt.plot(Urange,E_Gap_u,'-b',label='$Gap_u$')
plt.plot(Urange,E_Gap_d,'--r',label='$Gap_d$')
plt.xlabel('U(eV)')
plt.ylabel('Single-particle gap $\Delta_{sp}$')
plt.legend(loc='upper left')
#plt.grid(True)
plt.show()

plt.subplot(1,1,1)
plt.plot(Urange,Magnetization_a,'-b',label='$M_a$')
plt.plot(Urange,Magnetization_b,'--r',label='$M_b$')
plt.xlabel('U(eV)')
plt.ylabel('Magnetization')
plt.legend(loc='upper left')
#plt.grid(True)

#plt.tight_layout()
#plt.savefig( 'Energy_Gap_Ms_N'+str(N)+'.pdf')
#plt.savefig( 'KspaceResult.pdf')

plt.show()
