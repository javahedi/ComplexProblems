import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

plt.rcParams['axes.linewidth']     = 1.
plt.rcParams['xtick.major.size']   = 3
plt.rcParams['xtick.minor.size']   = 2
plt.rcParams['ytick.major.size']   = 3
plt.rcParams['ytick.minor.size']   = 2
plt.rcParams['xtick.bottom']       = True
plt.rcParams['xtick.top']          = True
plt.rcParams['ytick.left']         = True
plt.rcParams['ytick.right']        = True
plt.rcParams['xtick.direction']    = 'in'
plt.rcParams['ytick.direction']    = 'in'
plt.rcParams['font.family']        = 'sans-serif'
plt.rcParams['lines.linewidth']    = 2


#pathMSW  = 'MSW/chain/'
pathED   = 'ED/spin1'
pathFTLM = 'FTLM/spin1/chi'
pathKPM  = 'KPM/spin1/'
pathDMRG = 'DMRG/spin1/Roman_data2'
pathMSW  = 'MSW/chain'



colors = ['r','g','b']
styles = ['solid', 'dashed', 'dotted']
K_list = ['00','05','10']


fig, ax = plt.subplots(1,1,
            sharey=False,
            sharex=True,
            figsize=(5,3.5))



#######  MSM ##############
#Temperature,   <E>     F
#=================================================
msw00 = np.loadtxt(f'{pathMSW}/free_S1.0_K0.0_v2.txt')
msw05 = np.loadtxt(f'{pathMSW}/free_S1.0_K0.5_v2.txt')


#######  FTLM ##############
#Temperature,  F    <E>     <E^2>    Cv    entropy   chi    chi2   dE_LT      dE2_LT
#=================================================
ftlm00 = np.loadtxt(f'{pathFTLM}/L14D00/thermodynamic_FTLM_L14_K0.0.txt',skiprows=2)
ftlm05 = np.loadtxt(f'{pathFTLM}/L14D05/thermodynamic_FTLM_L14_K0.5.txt',skiprows=2)
ftlm10 = np.loadtxt(f'{pathFTLM}/L14D10/thermodynamic_FTLM_L14_K1.0.txt',skiprows=2)


#######  KPM ##############
#Temperature,  F    <E>     <E^2>    Cv    Chi    entropy
#=================================================
kpm00 = np.loadtxt(f'{pathKPM}/L10D00/thermodynamic_KPM_L10_K200.txt',skiprows=2)
kpm05 = np.loadtxt(f'{pathKPM}/L10D05/thermodynamic_KPM_L10_K200.txt',skiprows=2)
kpm10 = np.loadtxt(f'{pathKPM}/L10D10/thermodynamic_KPM_L10_K200.txt',skiprows=2)

#######  dmrg ##############
#T    β    c    e    chi    s
#=================================================
dmrg00 = np.loadtxt(f'{pathDMRG}/thermodyn_L=50_Kz=0_D=3_beta=100_SOFT=1.dat',skiprows=1)
dmrg05 = np.loadtxt(f'{pathDMRG}/thermodyn_L=50_Kz=-0.5_D=3_beta=100_SOFT=1.dat',skiprows=1)
dmrg10 = np.loadtxt(f'{pathDMRG}/thermodyn_L=50_Kz=-1_D=3_beta=100_SOFT=1.dat',skiprows=1)




ax.set_xlim(0,5)
ax.set_ylim(0,1.8)


ax.plot(msw00[:,0], -msw00[:,2]+1.4,  '-b', lw=1.0,label='MSW')
ax.plot(dmrg00[:,0], -dmrg00[:,3], label='DMRG,~L=40', lw=0.5,
                color='m', marker='o', markersize=4,
                markerfacecolor='none', markeredgecolor='m')
ax.plot(ftlm00[:,0], -(0.05+ftlm00[:,2]/14), label='FTLM,~~L=18', lw=0.5,
                color='c', marker='^', markersize=4,
                markerfacecolor='none', markeredgecolor='c')
ax.plot(kpm00[:,0], -(0.05+kpm00[:,2]/10), label='KPM,~~~~L=18', lw=0.5,
                color='y', marker='s', markersize=4,
                markerfacecolor='none', markeredgecolor='y')


ax.plot(msw05[:,0], -msw05[:,2]+0.35,  '-b', lw=1.0, )
ax.plot(dmrg05[:,0], -dmrg05[:,3], '-', lw=0.5,
                color='m', marker='o', markersize=4,
                 markerfacecolor='none', markeredgecolor='m')
ax.plot(ftlm05[:,0], -(0.05+ftlm05[:,2]/14), '-', lw=0.5,
                color='c', marker='^', markersize=4,
                 markerfacecolor='none', markeredgecolor='c')
ax.plot(kpm05[:,0], -(0.05+kpm05[:,2]/10),  lw=0.5,
                color='y', marker='s', markersize=4,
                markerfacecolor='none', markeredgecolor='y')




ax.minorticks_on()

ax.set_ylabel(r'$E(T)/N$',fontsize=14)
ax.set_xlabel(r'$T[1/J]$',fontsize=16)
ax.legend(ncol=1,loc='upper right',frameon=True,fontsize=10)
ax.text(0.38,0.2,r'$K=0.0$',transform=ax.transAxes,
            fontsize=14, rotation=0, color='k')
ax.text(0.38,0.6,r'$K=0.5$',transform=ax.transAxes,
            fontsize=14, rotation=0, color='k')

plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.99,
                    top=0.95,
                    wspace=0.0,
                    hspace=0.1)

plt.savefig('InternalEnergy.pdf', bbox_inches='tight')
plt.show()


