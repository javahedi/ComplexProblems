"""
Authors:   Javad + Guy
Date : June  2018  
Project:  The Hubbard model on Twisted Bilayer Graphene (TBG) in k_space
"""
import numpy as np
#import matplotlib.pyplot as plt
import time
import SelfConsistent
#import Magnetization2DPlot

##########################################
##########################################
##########################################
##########################################
#       Set parameters here
##########################################
##########################################
##########################################
##########################################
Dim = 1588    # number of atoms in twisted bilayer grahene
#U = 0.2        # Hubbard parameter in unit of (eV)
beta=1e6    # T=1e-3  temperature
Nx = 15     # number of k point in kx direction
Ny = 15     # number of k point in kx direction
NN = Nx *Ny  # number of unit cell
ne = 1       # number of electron per unit cell

################################################
#######         initial value           ########
################################################        
#site, Density_up, Density_dn = np.loadtxt('Initial_Density.dat', unpack=True)


#Density_up =   np.random.rand(Dim)
#Density_dn =   np.random.rand(Dim)

###  FM ####
Density_up = np.ones(Dim)
Density_dn = np.zeros(Dim)

#for i in range(Dim):
#    if np.mod(i,2)==0:
#       Density_up[i] = 1
#    else:
#       Density_dn[i] = 1




Mz_max = []
Mz_ave = []
E_GState = []
E_Gap_u = []
E_Gap_d = []

U  = np.loadtxt('Hubard_U.dat', unpack=True)
################################################
#######  Call self-consistent           ########
################################################
start = time.time()   
itr, Density_up, Density_dn, EG, gap_u, gap_d= \
	    SelfConsistent.SelfCon(Dim,U,Density_up,Density_dn,Nx,Ny,beta)

stop = time.time()

f5 = open('Dens_Nx'+str(Nx)+'_Ny'+str(Ny)+'_'+str(U)+'_FM.dat', 'w')

for i in range(Dim):
    f5.write("%i  %f   %f\n"  %(i+1, Density_up[i], Density_dn[i]))
	
f5.close()
################################################
#######      On-Site Magnetization   ###########
################################################    
Sz=np.subtract(Density_up,Density_dn)/2
Mz_ave.append(sum(abs(Sz))/Dim)
Mz_max.append(max(abs(Sz)))
################################################
#######      Single-Site Gap         ###########
################################################    
E_Gap_u.append(gap_u)
E_Gap_d.append(gap_d)
E_GState.append(EG)

f1 = open('Magnetization', 'w')
f2 = open('E_Gap_u.dat', 'w')
f3 = open('E_Gap_d.dat', 'w')
f4 = open('E_GState.dat', 'w')


f1.write("%f  %f   %f\n"  %(U, Mz_max[-1], Mz_ave[-1]))
f2.write("%f  %f\n"  %(U, E_Gap_u[-1]))
f3.write("%f  %f\n"  %(U, E_Gap_d[-1]))
f4.write("%f  %f\n"  %(U, E_GState[-1]))

print ('U = '+ str(U) +   ' ,itr = ' + str(itr) +\
   ', Mz_max = ' + str(Mz_max[-1]) +', Mz_ave = ' +str(Mz_ave[-1]) + ', Time = ' + str(stop - start))

f1.close()
f2.close()
f3.close()
f4.close()

