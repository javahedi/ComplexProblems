#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 14:09:29 2019

@author: javadvahedi
"""

import numpy as np
import os


a=[0.6,1.0,1.4,1.8,2.2,2.6]


os.chdir('job')
path0 = os.getcwd()

os.system('rm *.txt')

sxx_ave=[]
dist0_10=[]
for i in range(len(a)):
    os.chdir('alpha'+str(i)+'')    
    path1 = os.getcwd()
    
    ave=0
    dist=[]
    for j in range(1001):
        
        os.chdir('random'+str(j)+'')
        n,m,d,sxx0,sxx1=np.loadtxt('Sxx_L180_N_18_a_'+str(a[i])+'.txt',unpack=True)
        ave +=sxx1
        dist.append(sxx1[10])
        os.chdir(path1)
    sxx_ave.append(ave/1000)
    dist0_10.append(dist)
        
    os.chdir(path0)
    
    file1 = open('Sxx_ave_a_'+str(a[i])+'.txt', 'w')
    for l in range(len(sxx_ave[-1])):
        file1.write(" %f \n"  %(sxx_ave[-1][l]))
    
    
    file2 = open('dist0_10_a_'+str(a[i])+'.txt', 'w')
    for l in range(len(dist0_10[-1])):
        file2.write(" %f \n"  %(dist0_10[-1][l]))
    
    
    
file1.close()
file2.close()
     
    
    