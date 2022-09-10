
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import numpy as np
from numpy import pi, transpose, exp
import structure
import Hamiltonian
import Strain
import SelfConsistent
import Mayavi_New
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True


#%%
################################################
#######      Hamiltonian parameter      ########
################################################
t1 = -1.22  # eV
t2 =  3.665
t3 = -0.205
t4 = -0.105
t5 = -0.055
delta = 0.0
Hopping = [t1, t2, t3, t4, t5]
e0 = delta * abs(t1)
u = 1
U = u*abs(t1)
Nk = 6

f1 = open('Density_U_'+str(u)+'.txt', 'w')
f2 = open('EState_U_'+str(u)+'.txt', 'w')



#%%
################################################
#######      Strain parameter           ########
################################################
ex, ey, ez = 0.0, 0.0, 0.0

#%%
################################################
#######        Call  Hopping            ########
################################################
Hopping = Strain.strain(ex, ey, ez, Hopping)

#%%
################################################
#######      Lattice parameter          ########
################################################
Nx = 4   #Number of unit cell in x-direction 
Ny = 20  #Number of unit cell in y-direction 
pp = Ny//2

#%%
################################################
#######      Creating lattice           ########
################################################
Size = structure.struc(Nx,Ny)


#%%
################################################
#######          Call Hamiltonian       ########
################################################
Ham = Hamiltonian.Hamil(Hopping,e0,'Lattice_'+str(Size)+'.txt')

#%%
################################################
#######         Layer Hamiltonian       ########
################################################
Dim =  4 * Nx * pp   # each unit cell has 4 atom
H00 =  Ham[0:Dim,0:Dim]
H10 =  Ham[0:Dim,Dim:2*Dim]

#%%
################################################
#######    Initial configuration        ########
################################################
initial = 'FM'

if initial=='Random':
    Dens_up = np.random.rand(Dim)
    Dens_dn = np.random.rand(Dim)

if initial=='AF':
    Dens_up = np.empty(Dim)
    Dens_dn = np.empty(Dim)
    Dens_up[::2]  = 1
    Dens_up[1::2] = 0
    Dens_dn[::2]  = 0
    Dens_dn[1::2] = 1

if initial=='FM':
    Dens_up = np.ones(Dim)
    Dens_dn = np.zeros(Dim)

#%%
################################################
#######    self-consistent procedure    ########
################################################
Dens_up, Dens_dn, E_GState = \
SelfConsistent.SelfCon(Dim,U,Dens_up,Dens_dn,H00,H10,Nk,pp)

#%%
################################################
############    Magnetization           ########
################################################
Sz = (Dens_up-Dens_dn)/2
Mz_ave = sum(abs(Sz))/Dim
Mz_max = max(abs(Sz))


#%%
################################################
#######             Save data           ########
################################################
for i in range(Dim):
    f1.write("%i  %f  %f\n"  %(i+1, Dens_up[i], Dens_dn[i]))
f1.close()

f2.write("%f  %12.10f\n"  %(U, E_GState))
f2.close()

#%%
################################################
#######             Plot data           ########
################################################
plt.figure()
plt.plot(np.arange(Dim),Sz,'or')
plt.plot(np.arange(Dim),Sz,'-k')
plt.ylabel(r'$S^z$')
plt.xlabel('Site')
plt.title(r'$U/|t_1|=$' + str(u) + 'eV,  $E_G$='+str(E_GState)+'')
#plt.savefig('Sz.pdf')

Mayavi_New.Spatial(Nx,Ny,Size,Dim,u)




#%%
######################################
#######   Layer  Hamiltonian  ########
###################################### 

Dim = 4 * Nx  # each unit cell has 4 atom 
H00 =  Ham[0:Dim,0:Dim]
H10 =  Ham[0:Dim,Dim:2*Dim]

Huu = Hdd = H00
np.fill_diagonal(Huu,U*Dens_dn[:Dim])
np.fill_diagonal(Hdd,U*Dens_up[:Dim])


#%%
######################################
##########   Bandstructure    ########
######################################

Nk = 101
b = 3.32

EnergyBand_u=np.zeros((Nk,Dim))
EnergyBand_d=np.zeros((Nk,Dim))
ky=np.linspace(-pi/b,pi/b,Nk)
yy = np.zeros(Nk)

for k in range(Nk):

    HH = Huu + H10.T.conj() * exp(-1j*ky[k]*b)+\
               H10          * exp( 1j*ky[k]*b)
               
    en_u, vec_u = np.linalg.eig(HH)
    EnergyBand_u[k,:] = np.sort(np.real(en_u))

    HH = Hdd + H10.T.conj() * exp(-1j*ky[k]*b)+\
               H10          * exp( 1j*ky[k]*b)
               
    en_d, vec_d = np.linalg.eig(HH)    
    EnergyBand_d[k,:] = np.sort(np.real(en_d))
    
np.savetxt('Band_up_'+str(U)+'.txt', EnergyBand_u)
np.savetxt('Band_dn_'+str(U)+'.txt', EnergyBand_d)


#%% 
fig, ax1 = plt.subplots(figsize=(6,7))
for i in range(Dim):
    if min(abs(EnergyBand_u[:,i]))> 0.4:
        ax1.plot(ky,EnergyBand_u[:,i],'-b')
    else:
        ax1.plot(ky,EnergyBand_u[:,i],'-r')

ax1.plot(ky,yy,'--k')

rectangle = patches.Rectangle((-pi/b,0),         # (x,y)
                       2*pi/b,          # width
                       min(abs(EnergyBand_u[:,int(Dim/2)])),          # height
                       #fill=False,
                       edgecolor='None',
                       facecolor='yellow',
                       alpha = 1.0)
ax1.add_patch(rectangle)
plt.ylabel('$E-E_F[eV]$',fontsize = 20)
plt.xlabel('kb',fontsize = 20)
my_xticks = ['$-\\pi$', '0', '$\\pi$']
plt.xticks([-pi/b,0/b,pi/b], my_xticks, fontsize = 20)
plt.yticks(fontsize = 15)

plt.xlim(-pi/b,pi/b)
plt.ylim(-2,2)
#plt.title('$U= $'+str(U),fontsize = 20)


