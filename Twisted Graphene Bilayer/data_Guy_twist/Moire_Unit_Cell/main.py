#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 13:59:57 2018

@author: javad
"""
import numpy as np
import Configuration
import itertools
from numpy import pi, cos, sin, sqrt, log



def rotate_via_numpy(xy, radians):
    """Use numpy to build a rotation matrix and take the dot product."""
    x, y = xy
    c, s = np.cos(radians), np.sin(radians)
    j = np.matrix([[c, s], [-s, c]])
    m = np.dot(j, [x, y])

    return np.array([float(m.T[0]), float(m.T[1])])




a = 1      # c-c distance
d = 3 * a  # layers distance

r  = 1
m0 = 8

NN = 4 * (3*m0**2+3*m0*r+r**2/2)
N  = 4 * (3*m0**2+3*m0*r+r**2)
#theta = np.arccos(NN/N)*180/pi  in degree
theta = np.arccos(NN/N)


# fist's sheet vectors
a1 = a * np.array([ sqrt(3)/2, -1/2]) 
a2 = a * np.array([ sqrt(3)/2,  1/2])
d1 = (a1 + a2 )/3

# second's sheet vectors
aa1 = (cos(theta)-sin(theta)/sqrt(3)) * a1 + 2*sin(theta)/sqrt(3) * a2
aa2 = (cos(theta)+sin(theta)/sqrt(3)) * a2 - 2*sin(theta)/sqrt(3) * a1
dd1 = (aa1 + aa2 ) / 3


# superlattice vectors
R1 = m0*a1      + (m0+r)*a2
R2 = -(m0+r)*a1 + (2*m0+r)*a2

R1 = rotate_via_numpy(R1,0/6)
R2 = rotate_via_numpy(R2,0/6)

# the Moire period as a function of the twist angle
Lm = 1/np.sin(theta/2)/2
#Lm = np.linalg.norm(R1) or np.linalg.norm(R2)


#save data point
f1 = open('Lattice_r_'+str(r)+'m0_'+str(m0)+'.dat', 'w')


# generating lattice point with if conditoin to stauy in a big rectangel

site = 0
Nx , Ny =80, 80
NxList = np.arange(-Nx//2,(Nx+1)//2+1)
NyList = np.arange(-Ny//2,(Ny+1)//2+1)


xmax=20
ymax=16

# fist sheet data
for n,m in itertools.product(NxList, NyList):
    ra = n*a1 + m*a2
    rb = ra + d1
    if abs(ra[0])<xmax and abs(ra[1])<ymax:
        site += 1
        f1.write("%i  %f  %f  %f\n"  %(site, ra[0] , ra[-1], 0))   
        site += 1
        f1.write("%i  %f  %f  %f\n"  %(site, rb[0] , rb[-1], 0))

# second sheet data
for n,m in itertools.product(NxList, NyList):
    raa = n*aa1 + m*aa2 
    rbb = raa + dd1
    if abs(raa[0])<xmax and abs(raa[1])<ymax:
        site += 1
        f1.write("%i  %f  %f  %f\n"  %(site, raa[0] , raa[-1], d))    
        site += 1
        f1.write("%i  %f  %f  %f\n"  %(site, rbb[0] , rbb[-1], d))
        
f1.close()


theta = np.arccos(NN/N)*180/pi
Configuration.TwoDimPlot(r,m0,a,d,R1,R2,theta,N,xmax,ymax)
