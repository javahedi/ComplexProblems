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



pathKPM  = 'KPM/spin1/chi'
pathDMRG = 'DMRG/spin1/Roman_data2'




colors = ['r','g','b']
styles = ['solid', 'dashed', 'dotted']
K_list = ['00','05','10']


fig, ax = plt.subplots(2,1, sharey=False, sharex=True, figsize=(5,7))



ax0 = fig.add_axes([0.55, 0.60, 0.4, 0.2])
ax1 = fig.add_axes([0.55, 0.17, 0.4, 0.2])
ax0.minorticks_on()
ax1.minorticks_on()

#######  KPM ##############
#Temperature,  F    <E>     <E^2>    Cv    Chi    entropy
#=================================================

kpm05 = np.loadtxt(f'{pathKPM}/thermodynamic_KPM_L12_K0.5.txt',skiprows=2)



ax[0].set_xlim(0.01,1.1)
ax[1].set_xlim(0.01,1.1)
ax0.set_xlim(0.025,0.4)
ax1.set_xlim(0.025,0.4)
#ax[1].set_ylim(-0.01,1.1)
#ax[0].set_ylim(-0.001,0.2)


for j in range(100,900,100):
    kpmL8 = np.loadtxt(f'{pathKPM}/thermodynamic_KPM_L8_K0.5_NumberMoment{j}.txt',skiprows=2)
    
    ax[0].plot(kpmL8[:,0] , kpmL8[:,4]/8 ,label=r'$\mathcal{N}$='+str(j)+'', lw=1.5)
    ax0.loglog(kpmL8[:,0] , kpmL8[:,4]/8 , label=r'$\mathcal{N}$='+str(j)+'', lw=1.5)

    kpmL10 = np.loadtxt(f'{pathKPM}/thermodynamic_KPM_L10_K0.5_NumberMoment{j}.txt',skiprows=2)

    ax[1].plot(kpmL10[:,0] , kpmL10[:,4]/10 , label=r'$\mathcal{N}$='+str(j)+'', lw=1.5)
    ax1.loglog(kpmL10[:,0] , kpmL10[:,4]/10 , label=r'$\mathcal{N}$='+str(j)+'', lw=1.5)



ax[0].minorticks_on()
ax[1].minorticks_on()



ax[0].set_ylabel(r'$C_v(T)/L$',fontsize=16)
ax[1].set_ylabel(r'$C_v(T)/L$',fontsize=16)
ax[1].set_xlabel(r'$T[1/J]$',fontsize=16)
ax[1].legend(ncol=3,loc='upper right',frameon=True,fontsize=10)

ax[0].text(0.02,1.02,r'(a)',transform=ax[0].transAxes,
        fontsize=14, rotation=0, color='k')
ax[1].text(0.02,1.02,r'(b)',transform=ax[1].transAxes,
        fontsize=14, rotation=0, color='k')



ax0.set_ylabel(r'$C_v(T)/L$',fontsize=12)
ax1.set_ylabel(r'$C_v(T)/L$',fontsize=12)
ax0.set_xlabel(r'$T[1/J]$'  ,fontsize=10)
ax1.set_xlabel(r'$T[1/J]$'  ,fontsize=10)




ax[0].text(0.05,0.8,r'$L=14$',transform=ax[0].transAxes,
            fontsize=14, rotation=0, color='k')
ax[0].text(0.05,0.7,r'$K=0.5$',transform=ax[0].transAxes,
          fontsize=14, rotation=0, color='k')
ax[1].text(0.05,0.8,r'$L=18$',transform=ax[1].transAxes,
            fontsize=14, rotation=0, color='k')
ax[1].text(0.05,0.7,r'$K=0.5$',transform=ax[1].transAxes,
           fontsize=14, rotation=0, color='k')

plt.subplots_adjust(left=0.1, bottom=0.1,
                    right=0.99, top=0.95,
                    wspace=0.0, hspace=0.1)

plt.savefig('cv_compare_kmp.pdf', bbox_inches='tight')
plt.show()
