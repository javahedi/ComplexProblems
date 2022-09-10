
import scipy
import scipy.sparse as sparse
import numpy as np
#from numba import jit

#@jit(nopython=True)
class Hamiltonian(object):
    
    def __init__(self,op_list):
        self.op_list = op_list



    def gen_magnetization_squared(self):
        L = len(self.op_list)
        self.MM = self.op_list[0] * self.op_list[0]
        for i in range(1,L):
            self.MM = self.MM + self.op_list[i] * self.op_list[i]
            
        
    def gen_onsite_field(self,h_list):
        L = len(self.op_list)
        self.h = h_list[0] * self.op_list[0]
        for i in range(1,L): 
            self.h = self.h + h_list[i] * self.op_list[i]
        
    
    def gen_nn_int(self, J_list,bound,op_list2=[]):
        L = len(self.op_list)
    
        if op_list2 ==[]:
            op_list2 = self.op_list
            
        self.H = sparse.csr_matrix(self.op_list[0].shape)
        for ind,val in enumerate(bound):
            self.H = self.H +  J_list[ind] * self.op_list[val[0]] * op_list2[val[1]]
               
        
