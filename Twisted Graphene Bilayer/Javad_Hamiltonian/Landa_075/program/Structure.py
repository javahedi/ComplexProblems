#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 13:59:57 2018

@author: javad
"""
import numpy as np
#import Configuration
import itertools
from numpy import pi, cos, sin, sqrt, log


def UnitCell(m,n):

    #########################################    
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################

    a0 = 1.42 # c-c distance A
    a = 2.46  # Lattice constant a=np.sqrt(3) * a0
    d0 = 3.35  # Layers distance A



    #########################################    
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################

    def rotate(xy, radians):
        """Use numpy to build a rotation matrix and take the dot product."""
        x, y = xy
        c, s = cos(radians), sin(radians)
        j = np.matrix([[c, s], [-s, c]])
        m = np.dot(j, [x, y])

        return np.array([float(m.T[0]), float(m.T[1])])


    #########################################    
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################


    def graphene(m,n):      
        
        #  real primitive vector in unit of 
        a1 = a * np.array([ sqrt(3)/2, -0.5]) 
        a2 = a * np.array([ sqrt(3)/2,  0.5])
        d = (a1 + a2 )/3

        #  reciprocal primitive vector
        b1 = (2*pi/a) * (2/sqrt(3)) * np.array([   0.5, -sqrt(3)/2]) 
        b2 = (2*pi/a) * (2/sqrt(3)) * np.array([   0.5,  sqrt(3)/2])        
        
        # lattice vectors of the superlattice        
        L1 =  n*a1 +     m*a2
        L2 = -m*a1 + (m+n)*a2

        # lattice vectors of the reciprocal superlattice
        G1 = (   -m*b1 + n*b2 )/(n*n+m*m+m*n)
        G2 = ((m+n)*b1 + m*b2 )/(n*n+m*m+m*n)
        
        Nx , Ny = 70, 70
        NxList = np.arange(-Nx,(Nx+1))
        NyList = np.arange(-Ny,(Ny+1))
        
        
        x=[]
        y=[]
        site=[]
        cte=0
        for i,j in itertools.product(NxList, NyList):
            ra = i*a1 + j*a2        
            cte +=1
            x.append(ra[0])
            y.append(ra[1])
            site.append(cte)
            rb = ra + d
            cte +=1
            x.append(rb[0])
            y.append(rb[1])
            site.append(cte)
        


        return site,x,y,L1,L2, G1,G2

    #########################################    
    #########################################
    #########################################
    #########################################
    #########################################
    #########################################

    
    r1 = m**2 + n**2 + 4*m*n
    r2 = m**2 + n**2 + m*n
    #theta = np.arccos(NN/N)*180/pi  in degree
    theta = np.arccos(r1/r2/2)

    site,x,y,L1,L2, G1,G2 = graphene(m,n)

    ##################################
    #######   parallelogram  #########
    ##################################
    a0,b0 = L1,L2
    a1 = np.array([-a0[1],a0[0]]) 
    b1 = np.array([-b0[1],b0[0]])

    a2 = np.sign(np.dot(a1, b0)) * a1
    b2 = np.sign(np.dot(a0, b1)) * b1
    ##################################


    filename = 'unit_Cell_m_'+str(m)+'_n_'+str(n)+'.dat'
    f1 = open(filename, 'w')

    cte=1
    # layer one vectors
    for i in range(len(site)):    
        xy1 = rotate(np.array([x[i],y[i]]),0)+(L1+L2)/2
        con1= 0<np.dot(a2, xy1)<=np.dot(a2, b0)
        con2= 0<np.dot(b2, xy1)<=np.dot(a0, b2)
        if con1 & con2:
            f1.write("%i  %f  %f  %f\n"  %(cte, xy1[0] , xy1[1], 0))
            cte +=1
    
    # layer two vectors    
    for i in range(len(site)):
        xy2 = rotate(np.array([x[i],y[i]]),theta)+(L1+L2)/2
        con1= 0<=np.dot(a2, xy2)<=np.dot(a2, b0)
        con2= 0<np.dot(b2, xy2)<=np.dot(a0, b2)
        if con1 & con2:
            f1.write("%i  %f  %f  %f\n"  %(cte, xy2[0] , xy2[1], d0))
            cte +=1
        
    
    
    f1.close()


    theta = np.arccos(r1/r2/2)*180/pi    

    #Configuration.TwoDimPlot(n,m,1.42,d0,L1,L2,theta)
 

    return filename,L1,L2,G1,G2



