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



def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

colors = ['r','b','g']
styles = ['solid', 'dashed', 'dotted']
S_list = ['0.0','0.5']

'''
t05,l05 = np.loadtxt('lambda_S05.txt',unpack=True)
t1,l1  = np.loadtxt('lambda_S10.txt',unpack=True)
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
'''

################################
#  Free energy
################################
t00,l00 = np.loadtxt('lambda_S1.0_K0.0.txt',unpack=True)
t05,l05 = np.loadtxt('lambda_S1.0_K0.5.txt',unpack=True)

t00,f00 = np.loadtxt('free_S1.0_K0.0.txt',unpack=True)
t05,f05 = np.loadtxt('free_S1.0_K0.5.txt',unpack=True)

t00_,c00 = np.loadtxt('CV_S1.0_K0.0.txt',unpack=True)
t05_,c05 = np.loadtxt('CV_S1.0_K0.5.txt',unpack=True)


t00v2, eg00v2, f00v2 = np.loadtxt('free_S1.0_K0.0_v2.txt',unpack=True)
t05v2, eg05v2, f05v2 = np.loadtxt('free_S1.0_K0.5_v2.txt',unpack=True)



tem_list     = [t00,t05]
free_list    = [f00,f05]
free_list_v2 = [eg00v2,eg05v2]
tem_list_v2  = [t00v2,t05v2]

lamda_list = [l00,l05]
tem_list_  = [t00_,t05_]
c_list     = [c00,c05]


fig, ax = plt.subplots(figsize=(5,4))
ins     = fig.add_axes([0.22,0.5,0.34,0.36])
ax.minorticks_on()
ins.minorticks_on()
for id, K in enumerate([0.0,0.5]):   
    #ax.plot(tem_list[id],free_list[id],color=colors[id],ls=styles[id],lw=1,label=r'$K='+str(K)+'$')
    #ax.plot(tem_list_v2[id],free_list_v2[id],'.k')
    ax.plot(tem_list[id],lamda_list[id],color=colors[id],ls=styles[id],lw=1,label=r'$K='+str(K)+'$')
    

#########################
#    Spectrum
#########################
n = len(tem_list[0][0:250:10])

colors = plt.cm.jet(np.linspace(0,1,n))
cte =0
for T, l in zip(tem_list[0][0:250:10],lamda_list[0][0:250:10]):
    J=1.0
    K=0.0
    ###  mesh  over kx, ky
    nx = 501
    kx = np.linspace(-np.pi/2, np.pi/2, nx)
    N = nx
    gamma = np.cos(kx)
    p1 = (2*J + K)
    p0 = 2*J
    ek = np.sqrt((p1 * l)**2-(p0 * gamma)**2)
    ins.plot(kx/np.pi,ek+0.1,color=colors[cte],lw=1)
    cte +=1
    #print(T,l,abs(min(ek)-0.5))
#plt.grid('on')
#plt.show()


#ins.axhline(y=0.5, color='k', lw=1, linestyle='-')
#ax.axvline(x=0.6,  color='k', lw=1, linestyle='-')

#ins.yaxis.set_label_position("right")
#ins.yaxis.tick_right()
ins.set_ylabel(r'$\epsilon_k/J$',fontsize=12)
ins.set_xlabel(r'$k/\pi$',fontsize=12)
ins.set_xlim(-0.5,0.5)
ins.set_ylim(0.,4)

ax.set_ylabel(r'$\lambda[T]$',fontsize=16)
ax.set_xlabel(r'$T[1/J]$',fontsize=14)
ax.set_xlim(0.,3)
ax.set_ylim(0.7,1.8)
ax.legend(loc='lower right',frameon=False,fontsize=12)

ax.annotate("", xy=(1.74, 1.6), xytext=(1.74, 1.35),
            arrowprops=dict(arrowstyle="->",facecolor='k', lw=1))
ax.text(1.75,1.37,r"${\rm increasing}~~T$",rotation=90,fontsize=10)

#ins.text(0.18,0.6,r"$\Delta=0.5J$",rotation=0,fontsize=9)

plt.savefig('lambda.pdf', dpi=300, bbox_inches='tight')
plt.show()


######################################
######################################
######################################
######################################

S00 = -1.0 * np.diff(f00) 
S05 = -1.0 * np.diff(f05)

tem_list_S = [t00[:-1],t05[:-1]]
S_list     = [S00,S05]

C00 = np.diff(S00) * t00[:-2]
C05 = np.diff(S05) * t05[:-2]

tem_list_cv = [t00[:-2],t05[:-2]]
cv_list     = [C00,C05]

plt.figure(figsize=(5,4))
plt.minorticks_on()
plt.minorticks_on()
for id, K in enumerate([0.0,0.5]):
    plt.plot(tem_list_cv[id],1e2*cv_list[id],'-',label=r'K=' +str(K)+'')
    #plt.plot(tem_list_[id],   1.8*c_list[id],'--',label=r'K2='+str(K)+'')
    plt.plot(tem_list_cv[id], smooth(1e2*cv_list[id],9), '-', lw=2)
    
    columns = [ tem_list_cv[id], smooth(cv_list[id],9) ]
    np.savetxt(f'CV_S1.0_K{K}_smoot.txt',np.column_stack(columns))
    
    
plt.xlim(0.,3.)
#plt.ylim(-8.,0.)
plt.ylabel(r'$C_v/J$',fontsize=14)
plt.xlabel(r'$T/J$',fontsize=14)
plt.legend(loc='lower left',frameon=False,fontsize=12)
plt.savefig('Cv.pdf', dpi=300, bbox_inches='tight')
plt.show()





"""


################################
#  Entrophy
################################
t05 ,S05 = np.loadtxt('Entropy_S05.txt',unpack=True)
t1  ,S1  = np.loadtxt('Entropy_S10.txt',unpack=True)
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
"""
