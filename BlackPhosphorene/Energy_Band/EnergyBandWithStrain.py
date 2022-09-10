
from math import pi
import matplotlib.pyplot as plt
import numpy as np
import structure
import Hamiltonian
import Strain

import matplotlib as mpl
mpl.rcParams['text.usetex']        = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.family']       = 'sans'
#mpl.rcParams['font.family']        = 'DejaVu Sans'
mpl.rcParams['lines.linewidth']    = 2
mpl.rcParams['lines.markersize']   = 6
mpl.rcParams['font.size']          = 9 
mpl.rcParams['legend.fontsize']    = 14
mpl.rcParams['axes.titlesize']     = 14
mpl.rcParams['axes.labelsize']     = 14
mpl.rcParams['xtick.major.size']   = 3
mpl.rcParams['ytick.major.size']   = 3
mpl.rcParams['xtick.major.width']  = 1
mpl.rcParams['ytick.major.width']  = 1
mpl.rcParams['xtick.direction']    = 'in'
mpl.rcParams['ytick.direction']    = 'in'
mpl.rcParams['xtick.bottom']       = 'on'
mpl.rcParams['xtick.top']          = 'on'
mpl.rcParams['ytick.left']         = 'on'
mpl.rcParams['ytick.right']        = 'on'
mpl.rcParams['figure.titlesize']   = 14

# Hopping parameters
t1 = -1.22  # eV
t2 = 3.665
t3 = -0.205
t4 = -0.105
t5 = -0.055
delta = 0.0
Hopping = [t1, t2, t3, t4, t5]
e0 = delta*abs(t1)
U = 1
InitialValue = 0.5

# stran parameter
ex, ey = 0.0, 0.0
zplus = 1j*1e-6

Nx = 7  #Number of unit cell in x-direction 
Ny = 8  #Number of unit cell in y-direction 
Dim = structure.struc(Nx,Ny)




ezlist =  [0.0]
for ez in ezlist:
    Hopping = Strain.strain(ex, ey, ez, Hopping)
    ################################################
    #######  Call Hamiltonian  ########
    ################################################  
    Huu, Hdd, Density_up, Density_dn = \
             Hamiltonian.Hamil(Hopping,e0,InitialValue,U,'Lattice_'+str(Dim)+'.dat')

    N=len(Huu)
    ################################################
    #######    self-consistent procedure    ########
    ################################################


    beta = 1e5
    itr=0
    changeUp=1
    changeDn=1
    
    while ((changeUp > 1e-6) & (changeDn > 1e-6)):
        itr = itr +1
        Etot = []
        ################################################
        #######    Diagonalization           ###########
        ################################################
        eigs_u, vecs_u = np.linalg.eigh(Huu)
        eigs_d, vecs_d = np.linalg.eigh(Hdd)
        
        ################################################
        ####### Get chemical potential to    ###########
        ####### ensure n_up = N = n_down     ###########
        ################################################
       
        Etot.extend(eigs_u)
        Etot.extend(eigs_d)
        
        Etot.sort()
        fermi=(Etot[N]+Etot[N-1])/2
        
        ################################################
        ####### Finding New Density          ###########
        ################################################
        New_Density_up = np.zeros(N) 
        New_Density_dn = np.zeros(N)
          
        for j in range(N):
          prefac_u = 0.5*(1-np.tanh(beta*(eigs_u[j]-fermi)/2))
          prefac_d = 0.5*(1-np.tanh(beta*(eigs_d[j]-fermi)/2))
          New_Density_up[:] += prefac_u * vecs_u[:,j]*vecs_u[:,j]
          New_Density_dn[:] += prefac_d * vecs_d[:,j]*vecs_d[:,j]
          #for i in range(N):
          #      New_Density_up[i] += prefac_u * abs(vecs_u[i,j])**2#*vecs_u[i][j]
          #      New_Density_dn[i] += prefac_d * abs(vecs_d[i,j])**2#*vecs_d[i][j]
        
        
        p=0.7 # mixing parameter
        New_Density_up = p*Density_up+(1-p)*New_Density_up
        New_Density_dn = p*Density_dn+(1-p)*New_Density_dn
          
        changeUp=max(abs(New_Density_up-Density_up))
        changeDn=max(abs(New_Density_dn-Density_dn)) 
        
        Density_up = New_Density_up 
        Density_dn = New_Density_dn 
        
        ################################################
        ### Construction New Up and Down Hamiltonian  ##
        ################################################
        np.fill_diagonal(Huu, U*Density_dn)
        np.fill_diagonal(Hdd, U*Density_up)
        #for i in range(N):
        #    Huu[i][i]=U*Density_down[i] #-U*Density_down[i]*Density_up[i]*0.5
        #    Hdd[i][i]=U*Density_up[i] #-U*Density_down[i]*Density_up[i]*0.5
        
        print ('ez = ' + str(ez) + ', ' + ' Itr = ' + str(itr) +'' )
        ################################################
        #######      End of while loop     #############
        ################################################




    ######################################
    #######   Layer  Hamiltonian  ########
    ###################################### 
    NN = 4*Nx  # each unit cell has 4 atom
    H00_u = np.zeros((NN,NN))
    H10_u = np.zeros((NN,NN))
    H00_d = np.zeros((NN,NN))
    H10_d = np.zeros((NN,NN))
    nn = int(Ny/2)
    H00_u[0:NN,0:NN] =        Huu[nn*NN:(nn+1)*NN,nn*NN:(nn+1)*NN]
    H00_d[0:NN,0:NN] =        Hdd[nn*NN:(nn+1)*NN,nn*NN:(nn+1)*NN]

    H10_u[0:NN,0:NN] =        Huu[0:NN,NN:2*NN]
    H10_d[0:NN,0:NN] =        Hdd[0:NN,NN:2*NN]



    
    ######################################
    ##########   Bandstructure    ########
    ######################################
    
    Nk=201
    b = 3.32

    EnergyBand_u=np.zeros((Nk,NN))
    EnergyBand_d=np.zeros((Nk,NN))
    kx=np.linspace(-pi/b,pi/b,Nk)
    xx = np.zeros(Nk)

    for ii in range(Nk):

        HH_u=H00_u+np.transpose(H10_u).conj()*np.exp(-1j*kx[ii]*b)+\
               H10_u*np.exp(1j*kx[ii]*b)
        en_u, vec_u=np.linalg.eig(HH_u)    
        EnergyBand_u[ii,:]=np.sort(np.real(en_u))

        HH_d=H00_d+np.transpose(H10_d).conj()*np.exp(-1j*kx[ii]*b)+\
               H10_d*np.exp(1j*kx[ii]*b)
        en_d, vec_d=np.linalg.eig(HH_d)    
        EnergyBand_d[ii,:]=np.sort(np.real(en_d))
        
    np.savetxt('Band_up_ez_'+str(ez)+'.txt', EnergyBand_u)
"""/////////////////////////////////////////////////"""
"""/////////////////////////////////////////////////"""
 
#plt.figure(figsize=(4,4))
fig, ax1 = plt.subplots(figsize=(6,7))
for i in range(NN):
    if min(abs(EnergyBand_u[:,i]))> 0.4:
        ax1.plot(kx,EnergyBand_u[:,i],'-b')
    else:
        ax1.plot(kx,EnergyBand_u[:,i],'-r')

ax1.plot(kx,xx,'--k')

plt.ylabel('$E-E_F[eV]$',fontsize = 20)
plt.xlabel('kb',fontsize = 20)
my_xticks = ['$-\\pi$', '0', '$\\pi$']
plt.xticks([-pi/b,0/b,pi/b], my_xticks, fontsize = 20)
plt.yticks(fontsize = 15)

plt.xlim(-pi/b,pi/b)
plt.ylim(-3.3,3.0)
plt.title('$U= $'+str(U),fontsize = 20)




left, bottom, width, height = [0.7, 0.4, 0.18, 0.2]
ax2 = fig.add_axes([left, bottom, width, height])
for i in range(NN):
    if min(abs(EnergyBand_u[:,i]))> 0.4:
        ax2.plot(kx,EnergyBand_u[:,i],'-b')
    else:
        ax2.plot(kx,EnergyBand_u[:,i],'-r')

ax2.plot(kx,xx,'--k')
#plt.ylabel('$E-E_F[eV]$',fontsize = 20)
#plt.xlabel('kb',fontsize = 20)
my_xticks = ['$-\\pi$', '0', '$\\pi$']
plt.xticks([-pi/b,0/b,pi/b], my_xticks, fontsize = 10)
plt.yticks([])
plt.xlim(-pi/b,pi/b)
plt.ylim(-0.35,0.05)


plt.savefig('Band_Nx'+str(Nx)+'_Ny_'+str(Ny)+'_U_'+str(U)+'.png')
plt.show()






  



