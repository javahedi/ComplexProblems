import matplotlib.pyplot as pl
import numpy as np






EG = np.loadtxt('Ground_State_Energy.dat',  unpack=True)
Gap = np.loadtxt('Single_Particle_Gap.dat',  unpack=True)
Ma = np.loadtxt('Magnetization.dat',  unpack=True)



 

"""/////////////////////////////////////////////////"""
"""/////////////////////////////////////////////////"""

font = {'family': 'serif',
        #'color':  'darkred',
        'color':  'Black',
        'weight': 'bold',
        'size': 15,
        }

font2 = {'family': 'Times New Roman',
        #'color':  'darkred',
        'color':  'black',
        'weight': 'bold',
        'size': 15,
        }

"""/////////////////////////////////////////////////"""
"""/////////////////////////////////////////////////"""


Urange = np.arange(0,5,0.01)
Nt=np.size(Urange)
plot_list= [1, 2, 3]

f1 = open('Ground_State_Energy_new.dat', 'w')
f2 = open('Single_Particle_Gap_new.dat', 'w')
f3 = open('Magnetization_new.dat', 'w')
for i in range(Nt):
    f1.write("%f  %f\n"  %(Urange[i], EG[1,i] ))
    f2.write("%f  %f\n"  %(Urange[i], Gap[1,i] ))
    f3.write("%f  %f\n"  %(Urange[i], Ma[1,i] ))

f1.close()
f2.close()
f3.close()


fig1 = pl.figure(figsize=(5,5))
pl.subplot(3,1,1)
pl.plot(Urange,EG[1])
pl.xlabel('U/t', fontsize=15)
pl.ylabel('$E^G_{\\sigma}/N$', fontsize=10)

pl.subplot(3,1,2)
pl.plot(Urange,Gap[1])
pl.xlabel('U/t', fontsize=15)
pl.ylabel('Single particle gap', fontsize=10)
    
pl.subplot(3,1,3)
pl.plot(Urange,Ma[1])
pl.xlabel('U/t', fontsize=15)
pl.ylabel('Magnetization', fontsize=10)




pl.savefig( 'Graphene_Infinite.pdf')
pl.tight_layout() 
pl.show()
