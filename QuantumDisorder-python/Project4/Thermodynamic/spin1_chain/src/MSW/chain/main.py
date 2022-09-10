import sys
import numpy as np
import matplotlib.pyplot as plt

#############################
#   finding root with "bisection" method
#############################
def bisection(f, S, T, x_L, x_R, eps, return_x_list=False):
    """
    f   :  self-consistent function
    S   :  spin
    T   :  temperature
    x_L :  left  initial interval
    x_R :  right initial interval
    """
    f_L = f(x_L,S,T)
    if f_L*f(x_R,S,T) > 0:
        print("Error! Function does not have opposite \
                 signs at interval endpoints!")
        sys.exit(1)
    x_M = float(x_L + x_R)/2.0
    f_M = f(x_M,S,T)
    iteration_counter = 1
    if return_x_list:
        x_list = []

    while abs(f_M) > eps:
        if f_L*f_M > 0:   # i.e. same sign
            x_L = x_M
            f_L = f_M
        else:
            x_R = x_M
        x_M = float(x_L + x_R)/2
        f_M = f(x_M,S,T)
        iteration_counter += 1
        if return_x_list:
            x_list.append(x_M)
    if return_x_list:
        return x_list, iteration_counter
    else:
        return x_M, iteration_counter

#########################
#    self-consistent equation  Eq.6
#########################
def f(x,S,T):
    z = 2
    K = 0.5
    J = 1

    ###  mesh  over kx, ky
    nx = 51
    kx = np.linspace(-np.pi/2, np.pi/2, nx)
    N = nx
    gamma = np.cos(kx)
    p1 = (J * z + K) * S
    p0 = J * z * S
    ek = np.sqrt((p1 * x)**2-(p0 * gamma)**2)

    f0 = N * (2 * S + 1.0)
    coth = 1/np.tanh(ek/T/2.0)
    f1 = np.sum(coth*x*p1/ek)

    return f0-f1
##############



###  Define model paremeter #########
S = 1.0
T = 5 * S * (S+1)
# for D=0.0,0.5,1.0 ==> a=1.0, 0.8, 0.6
a = 0.8
b = 10
z = 2
K = 0.5
J = 1.0

f1 = open(f'lambda_S{S}_K{K}.txt','w')
landa = []
temp = []
for t in np.arange(0.01,T,0.02):
     
     solution, no_iterations = bisection(f, S, t, a, b, eps=1.0e-6)
     landa.append(solution)
     temp.append(t)
     f1.write('%f   %f\n'%(t,solution))
     f1.flush()
     #a = landa[-1]
f1.close()


#plt.figure()
#plt.plot(temp,landa,'-b')
#plt.show()





#########################
#    Free energy
#########################
f2 = open(f'free_S{S}_K{K}.txt','w')
free=[]
for T, l in zip(temp,landa):


    ###  mesh  over kx, ky
    nx = 51
    kx = np.linspace(-np.pi/2, np.pi/2, nx)
    N = nx

    gamma = np.cos(kx)
    p1 = (J * z + K) * S
    p0 = J * z * S
    ek = np.sqrt( (p1 * l)**2 - ( p0 * gamma)**2)

    E0 = -1.0 * S * (J+K) - S * (z*J+K) * (l-0.5)
    free.append(T*np.sum(np.log(2*np.sinh(ek/T/2.)))/N+E0)
    f2.write('%f   %f\n'%(T,-free[-1]))
    f2.flush()
f2.close()

plt.figure()
plt.plot(temp,free,'--b')
#plt.xlim([0.0,1.0])
#plt.ylim([0.0,1.0])
plt.show()

#########################
def bose(e,T):
    return 1.0/(np.exp(e/T)-1.0)
#########################



#########################
#    Internal energy
#########################
f4 = open(f'internal_S{S}_K{K}.txt','w')
internal=[]
for T, l in zip(temp,landa):


    ###  mesh  over kx, ky
    nx = 51
    kx = np.linspace(-np.pi/2, np.pi/2, nx)
    N = nx
    
    gamma = np.cos(kx)
    p1 = (J * z +K) * S
    p0 = J * z * S
    ek = np.sqrt((p1 * l)**2-(p0 * gamma)**2)

    E0 = 0.0#-1.0 * S * (J+K) - S * (z*J+K) * (l-0.5)
    internal.append(2.0*np.sum(bose(ek,T))/N+E0)
    f4.write('%f   %f\n'%(T,-internal[-1]))
    f4.flush()
f4.close()





#########################
#    Entropy
#########################

"""

f3 = open(f'entropy_S{S}_K{K}.txt','w')
entropy=[]
for T, l in zip(temp,landa):
    ###  mesh  over kx, ky
    nx = 51
    kx = np.linspace(-np.pi/2, np.pi/2, nx)
    N = nx

    gamma = np.cos(kx)
    p1 = (J * z + K) * S
    p0 = J * z * S
    ek = np.sqrt( (p1 * l)**2 - (p0 * gamma)**2 )

    nk1 = bose(ek,T) + 1.0
    nk  = bose(ek,T)
    entropy.append(np.sum(nk1*np.log(nk1)-nk*np.log(nk))/N)
    f3.write('%f   %f\n'%(T,entropy[-1]))
    f3.flush()
f3.close()


################################
#  CV
################################

cv = np.diff(np.array(entropy))*temp[:-1]
f5 = open(f'CV_S{S}_K{K}.txt','w')
cte=0
for T in temp[:-1]:
    f5.write('%f   %f\n'%(T,cv[cte]))
    f5.flush()
    cte+=1
f5.close()



################################
#  Chi
################################

f6 = open(f'chi_S{S}_K{K}.txt','w')
chi=[]
for T, l in zip(temp,landa):


    ###  mesh  over kx, ky
    nx = 51
    kx = np.linspace(-np.pi/2, np.pi/2, nx)
    N = nx

    gamma = np.cos(kx)
    p1 = (J * z + K) * S
    p0 = J * z * S
    ek = np.sqrt( (p1 * l)**2 - (p0 * gamma)**2 )

    
    chi.append(np.sum(bose(ek,T)*(bose(ek,T)+1.0))/N/3.0/T)
    f6.write('%f   %f\n'%(T,chi[-1]))
    f6.flush()
f6.close()


"""
