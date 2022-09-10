
import scipy
import scipy.sparse as sparse
import numpy as np



class Operator(object):

    def __init__(self,L):
        self.L = L

    
    def gen_op_total(self,op_list):
        #L = len(op_list)
        self.tot = op_list[0]
        for i in range(1,self.L):
            self.tot = self.tot + op_list[i]
    


    def gen_diagprojector(self,symvec,symval):
        ind = np.where(symvec==float(symval))
        dim = np.size(ind)
        self.P = sparse.lil_matrix((dim,len(symvec)))
        for j in range(dim):
            self.P[j,ind[0][j]] = 1.0
       
       
       
    def gen_s0sxsysz(self,ope='s0'):
        
        if ope=='s0':
            data = np.arange(2**self.L)
            row  = np.arange(2**self.L)
            col  = np.arange(2**self.L)
            s = sparse.csc_matrix((data,(row,col)),shape=(2**self.L,2**self.L))
        if ope=='sx':
            s = sparse.csc_matrix([[0., 1.],[1., 0.]])
        if ope=='sy':
            s = sparse.csc_matrix([[0.,-1j],[1j,0.]])
        if ope=='sz':
            s = sparse.csc_matrix([[1., 0],[0, -1.]])
        

        self.s_list = []
        for i_site in range(self.L):
            if i_site==0:
                S = s
            else:
                S= sparse.csc_matrix(np.eye(2))
            for j_site in range(1,self.L):
                if j_site==i_site:
                    S = sparse.kron(S,s, 'csc')
                else:
                    S = sparse.kron(S,np.eye(2),'csc')
            self.s_list.append(S)

    
