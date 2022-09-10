#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 14:31:25 2020

@author: javadvahedi
"""

import time
import resource

def check(start,date=False):
    print("="*50, flush=True)
    print(' elapset time         : '+str(time.time()-start)+' second', flush=True)
    checkmom = resource.getrusage(resource.RUSAGE_SELF)[2]/1e9
    print(' maximum memory usage : '+str(checkmom)+' GB', flush=True) 
    if date:
        print(' date                 :' ,time.strftime("%a,%d %b %Y %H:%M:%S",time.gmtime()), flush=True)
        print(" Run successful  :-)  ", flush=True)
        
    print("="*50, flush=True)
    return time.time()

