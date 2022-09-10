import Configuration

def struc(Nx,Ny):
   
    # lattice parameters
    r1 = 22.40  # nm
    r2 = 22.80
    r1x, r1y, r1z = 15.03, 16.60, 0.0
    r2x, r2y, r2z =  7.86,   0.0, 21.40
    a = 45.80
    b = 33.20
    
    #position of 4 atom in the nit cell
    atom_list = {0:[-r2x/2, 0,   0],\
                 1:[ r2x/2, 0, r2z],\
                 2:[-r2x/2-r1x, b/2,  0],\
                 3:[ r2x/2+r1x, b/2, r2z]}
    
    
    
    
    #Nx = 5    # number of unit cell in x direction
    #Ny = 2    # number of unit cell in y direction
    
    Dim = 4*Nx*Ny
    f1 = open('Lattice_'+str(Dim)+'.dat', 'w')
    
    site = 0
    for n in range(0,Ny):
        for m in range(0,Nx):
            for key in atom_list: # loop over 4 atom of unit cell 
                x =  m*a + atom_list[key][0]
                y =  n*b + atom_list[key][1]
                z =  atom_list[key][2]
                if y < (Ny * b)-20:
                    site += 1
                    f1.write("%i  %f  %f  %f  %i\n"  %(site, x , y, z, key))
               
            
    f1.close()
    
    #Configuration.TwoDimPlot(Dim,r1)
    return Dim
        
