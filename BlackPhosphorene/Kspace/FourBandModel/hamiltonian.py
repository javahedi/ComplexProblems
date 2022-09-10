from math import sqrt
import numpy as np

##############################
def Hamil(t,U,Density_up,Density_down,k):
    ax = 4.45       
    ay = 3.32
     

    ################################################
    #######  Call Up and Down  Hamiltonian  ########
    ################################################
    # elements of 4x4 up matrix
    f0 = 2*t[0]*np.exp(1j*ax*k[0]/(2*sqrt(3)))*np.cos(ay*k[1]/2)
    f1 = t[1]*np.exp(-1j*ax*k[0]/sqrt(3))
    f2 = 2*t[2]*np.exp(-1j*5*ax*k[0]/(2*sqrt(3)))*np.cos(ay*k[1]/2)
    f3 = 4*t[3]*np.cos(sqrt(3)*ax*k[0]/2)*np.cos(ay*k[1]/2)
    f4 = t[4]*np.exp(-1j*2*ax*k[0]/sqrt(3))
    
    g0 = f0.conj()
    g1 = f1.conj()
    g2 = f2.conj()
    g3 = f3.conj()
    g4 = f4.conj()
    
    
    Huu = np.array([[U*Density_down[0],             f0+f1,                f3,            f1+f4],
                    [            g0+g1, U*Density_down[1],                f1,               f3],
                    [               g3,                g1, U*Density_down[2],            f0+f2],
                    [            g1+g4,                g3,             g0+g2, U*Density_down[3]]])

    # elements of 4x4 down matrix
    Hdd = np.array([[U*Density_up[0],           f0+f1,              f3,          f1+f4],
                    [          g0+g1, U*Density_up[1],              f1,             f3],
                    [             g3,              g1, U*Density_up[2],          f0+f2],
                    [          g1+g4,              g3,           g0+g2, U*Density_up[3]]])

            
    ################################################N
    #######    Diagonalization           ###########
    ################################################
    eigs_u, vecs_u = np.linalg.eigh(Huu-(U/2)*np.eye(4))
    eigs_d, vecs_d = np.linalg.eigh(Hdd-(U/2)*np.eye(4))
    return eigs_u, vecs_u, eigs_d, vecs_d
