
import scipy
import scipy.sparse as sparse
import numpy as np


class Hamiltonian(object):
    
    def __init__(self,op_list):
        self.op_list = op_list


    def gen_onsite_field(self,h_list):
        L = len(self.op_list)
        self.h = h_list[0] * self.op_list[0]
        for i in range(1,L): 
            self.h = self.h + h_list[i] * self.op_list[i] 
        
    
    # generates \sum_i O_i O_{i+k} type interactions
    def gen_nn_int(self, k,J_list=[],op_list2=[],bc='obc'):
        L = len(self.op_list)
    
        if op_list2 ==[]:
            op_list2 = self.op_list
            
        self.H = sparse.csr_matrix(self.op_list[0].shape)
        if J_list == []:
            J_list = [1] * L
            
        Lmax = L if bc == 'pbc' else L-k
        for i in range(Lmax):
            self.H = self.H + J_list[i]*self.op_list[i]*op_list2[np.mod(i+k,L)]
                
        
