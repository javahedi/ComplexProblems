import numpy as np


def Hamil(t0,filename,m,n,L1,L2):
    #########################################    
    #########################################
    #########################################
    #########################################
    
   

    
    site, x, y, z = np.loadtxt(filename, unpack=True)
    #typ, x, y, z, site = np.loadtxt(filename, unpack=True)
    site = [ int(x) for x in site ]
    xL1,yL1,zL1=x+L1[0], y+L1[1], z
    xL2,yL2,zL2=x+L2[0], y+L2[1], z
    xL1L2,yL1L2,zL1L2=x+L1[0]+L2[0], y+L1[1]+L2[1], z
    xL1_L2,yL1_L2,zL1_L2=x+L1[0]-L2[0], y+L1[1]-L2[1], z


    N=len(site)       
         
           
    H00 = np.zeros((N,N)) 
       
    HL1   = np.zeros((N,N))
    HL2   = np.zeros((N,N))
    HL1L2 = np.zeros((N,N))
    HL1_L2= np.zeros((N,N))
        
    

        
    #########################################
    #########################################
    a0 = 1.418  # c-c distance Angestrum
    a  = 2.456  # Lattice constant a=np.sqrt(3) * a0
    d  = 3.349  # interlayer distance Angestrum
    #########################################
    #########################################
    
    #t0 = -2.7
    t1 = 0.48
    rc = 6.14
    lc = 0.265
    qs = 2.218*d
    qp = 2.218*a0
    
    
    
    for i in range(N):    
        for j in range(N):
            if i!=j:
                # interaction in each layer and between
                xij = x[i] - x[j]
                yij = y[i] - y[j]
                zij = z[i] - z[j]
                rij = np.sqrt(xij**2 + yij**2 + zij**2)
                n=(zij/rij)**2
                #if rij > a0/2:
                Fij = 1/(1+np.exp((rij-rc)/lc))        
                #if abs(t0*np.exp(qp*(1.0-rij/a0))*Fij) > 0.0001:
                if rij < 2.5*a:
                    H00[site[i]-1,site[j]-1] = (1.0-n)*t0*np.exp(qp*(1.0-rij/a0))*Fij+\
                                               n*t1*np.exp(qs*(1.0-rij/d))*Fij
	
        
        
                # interaction in L1 direction
                xij = x[i] - xL1[j]
                yij = y[i] - yL1[j]
                zij = z[i] - zL1[j]
                rij = np.sqrt(xij**2 + yij**2 + zij**2)
                
                n=(zij/rij)**2
                #if rij > a0/2:
                Fij = 1/(1+np.exp((rij-rc)/lc))            
                #if abs(t0*np.exp(qp*(1.0-rij/a0))*Fij) > 0.0001:
                if rij < 2.5*a:
                    HL1[site[i]-1,site[j]-1] = (1.0-n)*t0*np.exp(qp*(1.0-rij/a0))*Fij+\
                                               n*t1*np.exp(qs*(1.0-rij/d))*Fij


                # interaction in L2 direction
                xij = x[i] - xL2[j]
                yij = y[i] - yL2[j]
                zij = z[i] - zL2[j]
                rij = np.sqrt(xij**2 + yij**2 + zij**2)
                
                n=(zij/rij)**2
                #if rij > a0/2:
                Fij = 1/(1+np.exp((rij-rc)/lc))            
                #if abs(t0*np.exp(qp*(1.0-rij/a0))*Fij) > 0.0001:
                if rij < 2.5*a:
                    HL2[site[i]-1,site[j]-1] = (1.0-n)*t0*np.exp(qp*(1.0-rij/a0))*Fij+\
                                               n*t1*np.exp(qs*(1.0-rij/d))*Fij

               

                # interaction in L1 + L2 direction
                xij = x[i] - xL1L2[j]
                yij = y[i] - yL1L2[j]
                zij = z[i] - zL1L2[j]
                rij = np.sqrt(xij**2 + yij**2 + zij**2)
                
                n=(zij/rij)**2
                #if rij > a0/2:
                Fij = 1/(1+np.exp((rij-rc)/lc))            
                #if abs(t0*np.exp(qp*(1.0-rij/a0))*Fij) > 0.0001:
                if rij < 2.5*a:
                    HL1L2[site[i]-1,site[j]-1] = (1.0-n)*t0*np.exp(qp*(1.0-rij/a0))*Fij+\
                                                 n*t1*np.exp(qs*(1.0-rij/d))*Fij


               
                # interaction in L1 - L2 direction
                xij = x[i] - xL1_L2[j]
                yij = y[i] - yL1_L2[j]
                zij = z[i] - zL1_L2[j]
                rij = np.sqrt(xij**2 + yij**2 + zij**2)
                
                n=(zij/rij)**2
                #if rij > a0/2:
                Fij = 1/(1+np.exp((rij-rc)/lc))            
                #if abs(t0*np.exp(qp*(1.0-rij/a0))*Fij) > 0.0001:
                if rij < 2.5*a:
                    HL1_L2[site[i]-1,site[j]-1] = (1.0-n)*t0*np.exp(qp*(1.0-rij/a0))*Fij+\
                                                  n*t1*np.exp(qs*(1.0-rij/d))*Fij


               
    

    return H00,HL1,HL2,HL1L2,HL1_L2

