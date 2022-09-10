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


pathFTLM  = '../../FTLM/spin1/'
pathKPM  = '../../KPM/spin1/'

pathMy    = 'finiteT'
pathRoman = 'Roman_data2'



colors = ['r','g','b']
styles = ['solid', 'dashed', 'dotted']
K_list = ['00','05','10']


plt.figure(figsize=(5,4))



#######  KPM ##############
#Temperature,  F    <E> ...
#=================================================
kpm00 = np.loadtxt(f'{pathKPM}/L14D00/thermodynamic_KPM_L14_K0.0.txt',skiprows=2)
kpm05 = np.loadtxt(f'{pathKPM}/L14D05/thermodynamic_KPM_L14_K0.5.txt',skiprows=2)
kpm10 = np.loadtxt(f'{pathKPM}/L14D10/thermodynamic_KPM_L14_K1.0.txt',skiprows=2)

#######  FTLM ##############
#Temperature,  F    <E> ...
#=================================================
ftlm00 = np.loadtxt(f'{pathFTLM}/L14D00/thermodynamic_FTLM_L14_K0.0.txt',skiprows=2)
ftlm05 = np.loadtxt(f'{pathFTLM}/L14D05/thermodynamic_FTLM_L14_K0.5.txt',skiprows=2)
ftlm10 = np.loadtxt(f'{pathFTLM}/L14D10/thermodynamic_FTLM_L14_K1.0.txt',skiprows=2)


#######  My ##############
#Temperature,  beta    <E>
#=================================================
my05 = np.loadtxt(f'{pathMy}/K05/en_K0.5_cutoffE-10_Method_DensityMatrix.dat')
my10 = np.loadtxt(f'{pathMy}/K10/en_K1.0_cutoffE-10_Method_DensityMatrix.dat')




#######  Roman ##############
#T    β    c    e    chi    s
#=================================================
dmrg00 = np.loadtxt(f'{pathRoman}/thermodyn_L=50_Kz=0_D=3_beta=100_SOFT=1.dat',skiprows=1)
dmrg05 = np.loadtxt(f'{pathRoman}/thermodyn_L=50_Kz=-0.5_D=3_beta=100_SOFT=1.dat',skiprows=1)
dmrg10 = np.loadtxt(f'{pathRoman}/thermodyn_L=50_Kz=-1_D=3_beta=100_SOFT=1.dat',skiprows=1)


print(np.shape(my05))

plt.semilogx(my05[:,1],my05[:,2],'-b')
plt.semilogx(my10[:,1],my10[:,2],'-r')


plt.semilogx(dmrg00[:,0],dmrg00[:,3],'og')
plt.semilogx(dmrg05[:,0],dmrg05[:,3],'or')
plt.semilogx(dmrg10[:,0],dmrg10[:,3],'sb')

plt.semilogx(ftlm00[:,0],ftlm00[:,2]/14,'--k')
plt.semilogx(ftlm05[:,0],ftlm05[:,2]/14,'--k')
plt.semilogx(ftlm10[:,0],ftlm10[:,2]/14,'--k')


plt.semilogx(kpm00[:,0],kpm00[:,2]/14,'.g')
plt.semilogx(kpm05[:,0],kpm05[:,2]/14,'.r')
plt.semilogx(kpm10[:,0],kpm10[:,2]/14,'.b')




plt.ylim(-3,0)

#plt.subplots_adjust(left=0.1,
#                    bottom=0.1,
#                    right=0.99,
#                    top=0.95,
#                    wspace=0.0,
#                    hspace=0.1)
#plt.savefig('cv.pdf', bbox_inches='tight')
plt.show()


