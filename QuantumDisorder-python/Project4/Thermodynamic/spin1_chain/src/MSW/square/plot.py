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



t05,l05 = np.loadtxt('lambda_S05.txt',unpack=True)
t1,l1  = np.loadtxt('lambda_S1.txt',unpack=True)
t15,l15 = np.loadtxt('lambda_S15.txt',unpack=True)

tem_list   = [t05,t1,t15]
landa_list = [l05,l1,l15]

colors = ['r','g','b']
styles = ['solid', 'dashed', 'dotted']
S_list = ['1/2','1','3/2']


plt.figure(figsize=(5,4))
plt.minorticks_on()
plt.minorticks_on()
for id, S in enumerate([0.5,1.0,1.5]):
    
    plt.plot(tem_list[id]/(S*(S+1)),landa_list[id],color=colors[id],ls=styles[id],lw=2,label=r'S='+str(S_list[id])+'')

plt.xlim(0.,5.)
plt.ylim(0.9,2.5)
plt.ylabel(r'$\lambda(T)$',fontsize=16)
plt.xlabel(r'$T/JS(S+1)$',fontsize=14)
plt.legend(loc='upper left',frameon=False,fontsize=12)
plt.savefig('lambda.pdf', dpi=300, bbox_inches='tight')
plt.show()

################################
#  Free energy
################################
t05,f05 = np.loadtxt('FreeEnergy_S05.txt',unpack=True)
t1,f1  = np.loadtxt('FreeEnergy_S1.txt',unpack=True)
t15,f15 = np.loadtxt('FreeEnergy_S15.txt',unpack=True)

tem_list   = [t05,t1,t15]
free_list = [f05,f1,f15]


plt.figure(figsize=(5,4))
plt.minorticks_on()
plt.minorticks_on()
for id, S in enumerate([0.5,1.0,1.5]):
    
    plt.plot(tem_list[id]/(S*(S+1)),free_list[id]/(S*(S+1)),color=colors[id],ls=styles[id],lw=2,label=r'S='+str(S_list[id])+'')

plt.xlim(0.,5.)
plt.ylim(-8.,0.)
plt.ylabel(r'$F/JS(S+1)$',fontsize=14)
plt.xlabel(r'$T/JS(S+1)$',fontsize=14)
plt.legend(loc='lower left',frameon=False,fontsize=12)
plt.savefig('Free.pdf', dpi=300, bbox_inches='tight')
plt.show()




################################
#  Entrophy
################################
t05 ,S05 = np.loadtxt('Entropy_S05.txt',unpack=True)
t1  ,S1  = np.loadtxt('Entropy_S1.txt',unpack=True)
t15 ,S15 = np.loadtxt('Entropy_S15.txt',unpack=True)

tem_list = [t05,t1,t15]
ent_list = [S05,S1,S15]


plt.figure(figsize=(5,4))
plt.minorticks_on()
plt.minorticks_on()
for id, S in enumerate([0.5,1.0,1.5]):
    
    plt.plot(tem_list[id]/(S*(S+1)),ent_list[id],color=colors[id],ls=styles[id],lw=2,label=r'S='+str(S_list[id])+'')

plt.xlim(0.,5.)
plt.ylim(0.,1.7)
plt.ylabel(r'Entropy',fontsize=14)
plt.xlabel(r'$T/JS(S+1)$',fontsize=14)
plt.legend(loc='lower right',frameon=False,fontsize=12)
plt.savefig('Entropy.pdf', dpi=300, bbox_inches='tight')
plt.show()



################################
#  CV
################################

cv05 = np.diff(S05)*t05[:-1]
cv1  = np.diff(S1 )* t1[:-1]
cv15 = np.diff(S15)*t15[:-1]



tem_list = [t05[:-1],t1[:-1],t15[:-1]]
cv_list  = [cv05,cv1,cv15]

plt.figure(figsize=(5,4))
plt.minorticks_on()
plt.minorticks_on()
for id, S in enumerate([0.5,1.0,1.5]):
    
    plt.plot(tem_list[id]/(S*(S+1)),4*cv_list[id],color=colors[id],ls=styles[id],lw=2,label=r'S='+str(S_list[id])+'')

plt.xlim(0.,5.)
plt.ylim(0.,0.4)
plt.ylabel(r'Specific heat',fontsize=14)
plt.xlabel(r'$T/JS(S+1)$',fontsize=14)
plt.legend(loc='upper right',frameon=False,fontsize=12)
plt.savefig('SpeficiHeat.pdf', dpi=300, bbox_inches='tight')
plt.show()
