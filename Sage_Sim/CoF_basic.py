import copy
from sage.all import *


L = 2 # L transmitters
M = 2 # M relays
p = 3 # The prime number

Cores = 8 # The number of CPU cores used in parallel computing
DEBUG_H = False # When this value is True, the channel matrix H is set as certain matrices
P_Search_Alg = 'anneal' # 'brute', 'TNC', 'anneal'

def rate_computation(L, M, P_vec, alpha, H, A):
    r = zero_vector(RR, L)
    phi = [0]*M
    for i_m in range(0, M):
        sum_mis = 0
        for i_mis in range(0, L):
            sum_mis = sum_mis+(alpha[i_m]*H[i_m, i_mis]-A[i_m, i_mis])**2*P_vec[i_mis]
        phi[i_m] = (alpha[i_m])**2+sum_mis
    for i_l in range(0, L):
        if A.column(i_l).is_zero():
            r[i_l] = 0 # all coefficients are 0.
            raise Exception('It should be impossible')
        else:
            phi_max = 0
            for i_m in range(0, M):
                if A[i_m, i_l] != 0:
                    phi_max = max(phi[i_m], phi_max)
            r[i_l] = 0.5*log(max(1, P_vec[i_l]/phi_max), 2)
    return r

# Simple structure, for containing simulation result
class CoF_Sim_Result:
    def __init__(self, sum_rate_1, sum_rate_2):
        self.sum_rate_1 = sum_rate_1
        self.sum_rate_2 = sum_rate_2

class CoF_Sim_Result_for_Fixed_H_Fixed_P:
    def __init__(self, H_a, sum_rate_i_H, A, alpha_opt, P_vec):
        self.H_a = H_a
        self.sum_rate_i_H = sum_rate_i_H
        self.A = A
        self.alpha_opt = alpha_opt
        self.P_vec = P_vec

    def __repr__(self):
        return '---------------Fixed_H_Fixed_P---------------\n' \
            +'H_a = \n'+self.H_a.__str__()+'\n' \
            +'sum_rate_i_H = '+self.sum_rate_i_H.__str__()+'\n' \
            +'A = \n'+self.A.__str__()+'\n' \
            +'alpha_opt = '+self.alpha_opt.__str__()+'\n' \
            +'P_vec = '+self.P_vec.__str__()+'\n' \
            +'--------------------------------------------\n'
        
class CoF_Sim_Result_for_Fixed_H_Variable_P:
    def __init__(self, H_a, sum_rate_i_H_var, A, \
                 alpha_opt, P_vec_best_var):
        self.H_a = H_a
        self.sum_rate_i_H_var = sum_rate_i_H_var
        self.A = A
        self.alpha_opt = alpha_opt
        self.P_vec_best_var = P_vec_best_var

    def __repr__(self):
        return '-------------Fixed_H_Variable_P--------------\n' \
            +'H_a = \n'+self.H_a.__str__()+'\n' \
            +'sum_rate_i_H_var = '+self.sum_rate_i_H_var.__str__()+'\n' \
            +'A = \n'+ self.A.__str__()+'\n' \
            +'alpha_opt = '+self.alpha_opt.__str__()+'\n' \
            +'P_vec_best_var = '+self.P_vec_best_var.__str__()+'\n' \
            +'--------------------------------------------\n'