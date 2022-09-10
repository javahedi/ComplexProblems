"""
Authors:   Javad + Guy
Date : June  2018  
Project:  The Hubbard model on Twisted Bilayer Graphene (TBG) in k_space
"""
import numpy as np
import matplotlib.pyplot as plt
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
Dim = 868    # number of atoms in twisted bilayer grahene
#nk = 400     # number of k points
Nx=6
Ny=6

################################################
#######         initial value           ########
################################################        
#site, Dens= np.loadtxt('Initial_Density.dat', unpack=True)
#Density_up = 1 + 2 * Dens
#Density_dn = 1 - 2 * Dens

Density_up = np.zeros(Dim)
Density_dn = np.zeros(Dim)

for i in range(Dim):
    if np.mod(i,2)==0:
       Density_up[i] = 1
    else:
       Density_dn[i] = 1




U  = 0.8#np.loadtxt('Hubard_U.dat', unpack=True)
f1 = open('Magnetization_'+str(U)+'_.dat', 'w')
f2 = open('E_Gap_u_'+str(U)+'_.dat', 'w')
f3 = open('E_Gap_d_'+str(U)+'_.dat', 'w')
f4 = open('E_GState_'+str(U)+'_.dat', 'w')
f5 = open('Density_U_'+str(U)+'_.dat', 'w')

################################################
#######  Call self-consistent           ########
################################################
start = time.time()   
itr, Density_up, Density_dn, EG, gap_u, gap_d= \
	    SelfConsistent.SelfCon(Dim,U,Density_up,Density_dn,Nx,Ny)

stop = time.time()

for i in range(Dim):
    f5.write("%i  %f  %f\n"  %(i+1, Density_up[i], Density_dn[i]))
f5.close()

################################################
#######      On-Site Magnetization   ###########
################################################    
Sz=(Density_up-Density_dn)/2
Mz_ave = sum(abs(Sz))/Dim
Mz_max = max(abs(Sz))
Mz_1 = sum(Sz[0:Dim//2])/(Dim//2)
Mz_2 = sum(Sz[Dim//2:Dim])/(Dim//2)
Mz_1_2 = sum(Sz)/Dim
################################################
#######      Single-Site Gap         ###########
################################################    


f1.write("%f  %f   %f  %f   %f   %f\n"  %(U, Mz_max, Mz_ave,\
                                          Mz_1, Mz_2, Mz_1_2))
f2.write("%f  %f\n"  %(U, gap_u))
f3.write("%f  %f\n"  %(U, gap_d))
f4.write("%f  %f\n"  %(U, EG))



print ('U = '+ str(U) +   ' ,itr = ' + str(itr) +\
   ', Mz_max = ' + str(Mz_max) +', Mz_ave = ' +str(Mz_ave) + \
       ', Time = ' + str(stop - start))

f1.close()
f2.close()
f3.close()
f4.close()

