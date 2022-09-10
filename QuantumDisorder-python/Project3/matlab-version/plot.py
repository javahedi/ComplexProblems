import numpy as np
import matplotlib.pyplot as plt

alpha=1.0

Con=np.loadtxt('Con.txt',unpack=True)
SySy=np.loadtxt('SySy.txt',unpack=True)


fig = plt.figure(figsize=(4,4))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],facecolor='w')
ax.tick_params(direction='in',axis='both', color='k',labelsize=14,
               left='on', top='on', right='on', bottom='on')

ax.loglog(abs(SySy), "o-b")
x=np.arange(1,len(SySy))
ax.loglog(x,np.exp(-x),'-r')
ax.loglog(x,1/(x**alpha),'-c')
ax.axes.set_xlim([0.9,11])
ax.axes.set_xlabel(r"$|n-m|$",fontsize=18)
ax.axes.set_ylabel(r"$\langle S^y_nS^y_m\rangle$",fontsize=18)
plt.show()
