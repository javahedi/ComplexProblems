import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

from matplotlib import rc
import matplotlib.image as mpimg
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


pathDMRG1 = 'DMRG/spin1/Roman_data1'
pathDMRG2 = 'DMRG/spin1/Roman_data2'
pathDMRG3 = 'DMRG/spin1/finiteT'




colors = ['r','g','b']
styles = ['solid', 'dashed', 'dotted']
K_list = ['00','05','10']


fig, ax = plt.subplots(1,1, figsize=(5,3.5))



#ax0 = fig.add_axes([0.65, 0.62, 0.31, 0.16])
#ax1 = fig.add_axes([0.65, 0.22, 0.31, 0.16])
#ax0.minorticks_on()
#ax1.minorticks_on()

#######  dmrg ##############
#T    β    c    e    chi    s
#=================================================
dmrg1_00 = np.loadtxt(f'{pathDMRG1}/thermodyn_L=50_D=3_beta=100.dat',skiprows=1)
dmrg1_05 = np.loadtxt(f'{pathDMRG1}/thermodyn_L=50_Kz=-0.5_D=3_beta=100.dat',skiprows=1)
dmrg1_10 = np.loadtxt(f'{pathDMRG1}/thermodyn_L=50_Kz=-1_D=3_beta=100.dat',skiprows=1)


dmrg2_00 = np.loadtxt(f'{pathDMRG2}/thermodyn_L=50_Kz=0_D=3_beta=100_SOFT=1.dat',skiprows=1)
dmrg2_05 = np.loadtxt(f'{pathDMRG2}/thermodyn_L=50_Kz=-0.5_D=3_beta=100_SOFT=1.dat',skiprows=1)
dmrg2_10 = np.loadtxt(f'{pathDMRG2}/thermodyn_L=50_Kz=-1_D=3_beta=100_SOFT=1.dat',skiprows=1)



dmrg3_05 = np.loadtxt(f'{pathDMRG3}/K05/sus_K0.5_cutoffE-10_Method_DensityMatrix.dat')



ax.set_xlim(0.01,100)
#ax.set_ylim(0.0,3)



ax.semilogx(dmrg1_00[:,0] , dmrg1_00[:,4] , '-b',label=r'$K=0.0$, OBC', lw=1.5)
ax.semilogx(dmrg1_05[:,0] , dmrg1_05[:,4] , '-g',label=r'$K=0.5$, OBC', lw=1.5)
ax.semilogx(dmrg1_10[:,0] , dmrg1_10[:,4] , '-r',label=r'$K=1.0$, OBC', lw=1.5)


ax.semilogx(dmrg2_00[:,0] , dmrg2_00[:,4] , '.b', label=r'$K=0.0$, SBC')
ax.semilogx(dmrg2_05[:,0] , dmrg2_05[:,4] , '.g', label=r'$K=0.5$, SBC')
ax.semilogx(dmrg2_10[:,0] , dmrg2_10[:,4] , '.r', label=r'$K=1.0$, SBC')

ax.legend(ncol=1,loc='upper left',frameon=True,fontsize=8)



ax.minorticks_on()



ax.set_ylabel(r'$\chi(T)/L$',fontsize=14)
ax.set_xlabel(r'$T[1/J]$',fontsize=16)




#ax0.set_ylabel(r'$\chi(T)/N$',fontsize=12)
#ax0.set_xlabel(r'$T[1/J]$',fontsize=10)


ax.text(0.38,0.9,r'OBC',transform=ax.transAxes,
        fontsize=12, rotation=0, color='k')
ax.text(0.38,0.73,r'SBC',transform=ax.transAxes,
        fontsize=12, rotation=0, color='k')

'''
ax.text(0.2,0.5,r'OBC',transform=ax.transAxes,
            fontsize=12, rotation=0, color='k')
            
ax.annotate("", xy=(0.028, 12), xytext=(0.1, 16),
             arrowprops=dict(arrowstyle="->",facecolor='red', lw=1))
ax.annotate("", xy=(0.027, 6), xytext=(0.1, 16),
            arrowprops=dict(arrowstyle="->",facecolor='green', lw=1))
ax.annotate("", xy=(0.028, 1), xytext=(0.1, 16),
            arrowprops=dict(arrowstyle="->",facecolor='blue', lw=1))
'''

#ax.text(0.1,0.8,r'$L=40$',transform=ax.transAxes,
#        fontsize=14, rotation=0, color='k')


ax0 = fig.add_axes([0.55, 0.275, 0.4, 0.35])
ax0.minorticks_on()
ax0.set_ylabel(r'$\chi(T)/N$',fontsize=10)
ax0.set_xlabel(r'$T[1/J]$',fontsize=12)
ax0.set_ylim(0.0,2)
ax0.set_xlim(0.01,1.1)


ax0.semilogx(dmrg1_00[:,0] , dmrg1_00[:,4] , '-b',label=r'K=0.0', lw=1.5)
ax0.semilogx(dmrg1_05[:,0] , dmrg1_05[:,4] , '-g',label=r'K=0.5', lw=1.5)
ax0.semilogx(dmrg1_10[:,0] , dmrg1_10[:,4] , '-r',label=r'K=1.0', lw=1.5)


ax0.semilogx(dmrg2_00[:,0] , dmrg2_00[:,4] , '.b', label=r'K=0.0')
ax0.semilogx(dmrg2_05[:,0] , dmrg2_05[:,4] , '.g', label=r'K=0.5')
ax0.semilogx(dmrg2_10[:,0] , dmrg2_10[:,4] , '.r', label=r'K=1.0')



im = mpimg.imread('paper/Figure/Cartoon.png', format='png')
imagebox = OffsetImage(im, zoom=0.2)

ab = AnnotationBbox(imagebox, (8, 28), frameon=False)
ax.add_artist(ab)

#axin.axis('off')
#plt.draw()





plt.subplots_adjust(left=0.1, bottom=0.1,
                    right=0.99, top=0.95,
                    wspace=0.0, hspace=0.1)

plt.savefig('chi_compare_dmrg.pdf', bbox_inches='tight')
plt.show()


