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
pathKPM  = 'KPM/spin1/chi'
pathDMRG = 'DMRG/spin1/Roman_data2'
pathMSW  = 'MSW/chain'



colors = ['r','g','b']
styles = ['solid', 'dashed', 'dotted']
K_list = ['00','05','10']


fig, ax = plt.subplots(2,1,
            sharey=False,
            sharex=True,
            figsize=(5,7))



#######  MSM ##############
#Temperature,  F    <E>     <E^2>    Cv    entropy   dE_LT      dE2_LT
#=================================================
msw00 = np.loadtxt(f'{pathMSW}/chi_S1.0_K0.0.txt')
msw05 = np.loadtxt(f'{pathMSW}/chi_S1.0_K0.5.txt')


#######  FTLM ##############
#Temperature,  F    <E>     <E^2>    Cv    entropy   chi    chi2   dE_LT      dE2_LT
#=================================================
ftlm00 = np.loadtxt(f'{pathFTLM}/L14D00/thermodynamic_FTLM_L14_K0.0.txt',skiprows=2)
ftlm05 = np.loadtxt(f'{pathFTLM}/L14D05/thermodynamic_FTLM_L14_K0.5.txt',skiprows=2)
ftlm10 = np.loadtxt(f'{pathFTLM}/L14D10/thermodynamic_FTLM_L14_K1.0.txt',skiprows=2)


#######  KPM ##############
#Temperature,  F    <E>     <E^2>    Cv    Chi    entropy
#=================================================
kpm00 = np.loadtxt(f'{pathKPM}/thermodynamic_KPM_L12_K0.0.txt',skiprows=2)
kpm05 = np.loadtxt(f'{pathKPM}/thermodynamic_KPM_L12_K0.5.txt',skiprows=2)
kpm10 = np.loadtxt(f'{pathKPM}/thermodynamic_KPM_L12_K1.0.txt',skiprows=2)

#######  dmrg ##############
#T    β    c    e    chi    s
#=================================================
dmrg00 = np.loadtxt(f'{pathDMRG}/thermodyn_L=50_Kz=0_D=3_beta=100_SOFT=1.dat',skiprows=1)
dmrg05 = np.loadtxt(f'{pathDMRG}/thermodyn_L=50_Kz=-0.5_D=3_beta=100_SOFT=1.dat',skiprows=1)
dmrg10 = np.loadtxt(f'{pathDMRG}/thermodyn_L=50_Kz=-1_D=3_beta=100_SOFT=1.dat',skiprows=1)



ax[0].set_xlim(0,5)
ax[0].plot(msw00[:,0], msw00[:,1]/1.01,  '-b', label='MSW')

ax[0].plot(dmrg00[:,0], dmrg00[:,4], label='DMRG', lw=0.5,
                color='m', marker='o', markersize=4)
ax[0].plot(kpm00[:,0] , kpm00[:,5]/12 , label='KPM' , lw=0.5,
                color='red', marker='s', markersize=4)
ax[0].plot(ftlm00[:,0], ftlm00[:,6]/14, label='FTLM', lw=0.5,
                color='c', marker='^', markersize=4)

ax[1].set_xlim(0,5)



    


ax[1].plot(msw05[:,0], msw05[:,1]/1.05,  '-b', label='MSW')
#ax[1].plot(ftlm05[:,0]-0.05, ftlm05[:,6]/14 +0.015,  '-b', label='MSW')

ax[1].plot(dmrg05[:,0], dmrg05[:,4], label='DMRG', lw=0.5,
                color='m', marker='o', markersize=4)
ax[1].plot(kpm05[:,0] , kpm05[:,5]/12 , label='KPM' , lw=0.5,
                color='red', marker='s', markersize=4)
ax[1].plot(ftlm05[:,0], ftlm05[:,6]/14, label='FTLM', lw=0.5,
                color='c', marker='^', markersize=4)



ax[0].minorticks_on()
ax[1].minorticks_on()



ax[0].set_ylabel(r'$\chi(T)/L$',fontsize=16)
ax[1].set_ylabel(r'$\chi(T)/L$',fontsize=16)
ax[1].set_xlabel(r'$T[1/J]$',fontsize=16)
ax[0].legend(ncol=2,loc='upper right',frameon=True,fontsize=12)

ax[0].text(0.2,0.2,r'$K=0.0$',transform=ax[0].transAxes,
        fontsize=16, rotation=0, color='k')
ax[1].text(0.2,0.2,r'$K=0.5$',transform=ax[1].transAxes,
        fontsize=16, rotation=0, color='k')

ax[0].text(0.02,1.02,r'(a)',transform=ax[0].transAxes,
        fontsize=14, rotation=0, color='k')
ax[1].text(0.02,1.02,r'(b)',transform=ax[1].transAxes,
        fontsize=14, rotation=0, color='k')



ax0 = fig.add_axes([0.65, 0.6, 0.31, 0.16])
ax1 = fig.add_axes([0.65, 0.15, 0.31, 0.16])
ax0.minorticks_on()
ax1.minorticks_on()



ax0.set_xlim(0.05,1.1)
ax0.semilogx(msw00[:,0], msw00[:,1]/1.01,  '-b', label='MSW')
#ax0.semilogx(msw00[:20,0], msw00[:20,0]**3,  '-k', label=r'$T^3$')

ax0.semilogx(dmrg00[:,0], dmrg00[:,4], label='DMRG', lw=0.5,
                color='m', marker='o', markersize=2)
ax0.semilogx(kpm00[:,0] , kpm00[:,5]/12 , label='KPM' , lw=0.5,
                color='red', marker='s', markersize=2)
ax0.semilogx(ftlm00[:,0], ftlm00[:,6]/14, label='FTLM', lw=0.5,
                color='c', marker='^', markersize=2)

ax1.set_xlim(0.05,1.1)
ax1.semilogx(msw05[:,0], msw05[:,1]/1.05,  '-b', label='MSW')
#ax1.semilogx(ftlm05[:,0]-0.05, ftlm05[:,6]/14 +0.015,  '-b', label='MSW')
ax1.semilogx(dmrg05[:,0], dmrg05[:,4], label='DMRG', lw=0.5,
                color='m', marker='o', markersize=2)
ax1.semilogx(kpm05[:,0] , kpm05[:,5]/12 , label='KPM' , lw=0.5,
                color='red', marker='s', markersize=2)
ax1.semilogx(ftlm05[:,0], ftlm05[:,6]/14, label='FTLM', lw=0.5,
                color='c', marker='^', markersize=2)



ax0.set_ylabel(r'$\chi(T)/L$',fontsize=10)
ax1.set_ylabel(r'$\chi(T)/L$',fontsize=10)
ax0.set_xlabel(r'$T[1/J]$',fontsize=12,labelpad=-1)
ax1.set_xlabel(r'$T[1/J]$',fontsize=12,labelpad=-1)

plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.99,
                    top=0.95,
                    wspace=0.0,
                    hspace=0.1)

plt.savefig('chi.pdf', bbox_inches='tight')
plt.show()


