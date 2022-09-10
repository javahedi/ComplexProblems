import sys
import numpy as np 


#############################
#   finding root with "bisection" method
#############################
def bisection(f, S, T, x_L, x_R, eps, return_x_list=False):
    """
    f  : Eq.6
    S  : spin
    T  : temperature
    x_L : left initial interval
    x_R : right initial interval
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
    J = 1
    z = 4
    
    ###  mesh  over kx, ky
    nx, ny = 51, 51
    kx = np.linspace(-np.pi, np.pi, nx)
    ky = np.linspace(-np.pi, np.pi, ny)
    kx, ky = np.meshgrid(kx, ky)
    N = nx * ny

    gamma = np.cos(kx/2)*np.cos(ky/2)
    ek = J*z*S*np.sqrt(x**2-gamma**2)
    
    f0 = 2*N*(S+0.5)
    
    coth = 1/np.tanh(ek/T/2.0)
    f1 = np.sum(coth*x/np.sqrt(x**2-gamma**2))
    
    return f0-f1
#################################



###  Define model paremeter #########
S = 0.5
T = 5 * S * (S+1)
a = 1
b = 10

f1 = open('lambda_S05.txt','w')
landa = []
temp = []
for t in np.arange(0.1,T,0.05):
     solution, no_iterations = bisection(f, S, t, a, b, eps=1.0e-6)
     landa.append(solution)
     temp.append(t)
     f1.write('%f   %f\n'%(t,solution))
     f1.flush()
f1.close()   




#########################
#    Free energy  Eq.5
#########################
f2 = open('FreeEnergy_S05.txt','w')
free=[]
for T, l in zip(temp,landa):
    J = 1
    z = 4
    
    ###  mesh  over kx, ky
    nx, ny = 51, 51
    kx = np.linspace(-np.pi, np.pi, nx)
    ky = np.linspace(-np.pi, np.pi, ny)
    kx, ky = np.meshgrid(kx, ky)
    N = nx * ny

    gamma = np.cos(kx/2)*np.cos(ky/2)
    ek = J*z*S*np.sqrt(l**2-gamma**2)
    
    free.append(T*np.sum(np.log(2*np.sinh(ek/T/2.)))/N-J*z*l*S*(1+2*S)/2 + J*z*S**2/2)
    f2.write('%f   %f\n'%(T,free[-1]))
    f2.flush()
f2.close()




#########################
#    Entropy  Eq.5
#########################
def bose(e,T):
    return 1.0/(np.exp(e/T)-1.0)



f3 = open('Entropy_S05.txt','w')
entropy=[]
for T, l in zip(temp,landa):
    J = 1
    z = 4
    
    ###  mesh  over kx, ky
    nx, ny = 51, 51
    kx = np.linspace(-np.pi, np.pi, nx)
    ky = np.linspace(-np.pi, np.pi, ny)
    kx, ky = np.meshgrid(kx, ky)
    N = nx * ny

    gamma = np.cos(kx/2)*np.cos(ky/2)
    ek = J*z*S*np.sqrt(l**2-gamma**2)
    
    nk1 = bose(ek,T) + 1.0
    nk  = bose(ek,T)
    entropy.append(np.sum(nk1*np.log(nk1)-nk*np.log(nk))/N)
    f3.write('%f   %f\n'%(T,entropy[-1]))
    f3.flush()
f3.close()
    

