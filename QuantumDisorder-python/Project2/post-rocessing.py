import numpy as np
from scipy.stats import gmean
import os

n_random = 1000
L = 120
alpha = [1.0]

path0 = os.getcwd()
os.chdir('job')
path1 = os.getcwd()

for a in alpha:
    print( a, end='\r',flush=True)

    os.chdir('alpha'+str(a)+'')
    
    svl_ave, svs_ave = 0,0
    sxx_ave = np.zeros(L)
    typical = [[] for i in range(L)]
    for w in range(1,n_random):
        l, svl = np.loadtxt('Svn_RDM_length_'+str(w)+'.txt',unpack=True)
        svl_ave +=svl
        
        svs = np.loadtxt('Svn_RDM_site_'+str(w)+'.txt')
        svs_ave +=svs

        m, n,dist, sxx, sxx1 = np.loadtxt('Sxx_'+str(w)+'.txt',unpack=True)
        indx = dist.argsort()
        dist = dist[indx]
        sxx1 = sxx1[indx]
       
        
        for i in range(len(dist)):
            t = int(dist[i])
            sxx_ave[t] +=sxx1[i]
            typical[t].append(sxx1[i])
           
    sxx_typical = np.zeros(L)
    for id, val in enumerate(typical):    
        val = list(filter(lambda num: num != 0, val))
        sxx_typical[id] = gmean(val)

    sxx_typical = list(filter(lambda x: str(x) != 'nan', sxx_typical))
    
    svl_ave_sym = (svl_ave + svl_ave[::-1])/2
    os.chdir(path0)

    np.savetxt('svn_alpha_'+str(a)+'_length_'+str(L)+'_L.txt',svl_ave/n_random)
    np.savetxt('svn_alpha_'+str(a)+'_length_'+str(L)+'_L_sym.txt',svl_ave_sym/n_random)
    np.savetxt('svn_alpha_'+str(a)+'_length_'+str(L)+'_S.txt',svs_ave/n_random)
    np.savetxt('sxx_alpha_'+str(a)+'_length_'+str(L)+'.txt',sxx_ave/n_random)
    np.savetxt('typical_sxx_alpha_'+str(a)+'_length_'+str(L)+'.txt',sxx_typical)
    os.chdir(path1)
    
    


