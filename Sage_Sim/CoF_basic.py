import copy
from sage.all import *
import math

L = 2 # L transmitters
M = 2 # M relays
p = 3 # The prime number

Cores = 8 # The number of CPU cores used in parallel computing
DEBUG_H = False # When this value is True, the channel matrix H is set as certain matrices
P_Search_Alg = 'brute' # 'brute', 'TNC', 'anneal'
brute_number = 50
is_alternate = False # True or False
iter_H = 32
batch_H = 4

def alpha_find(h, P_mat, a):
    alpha_opt = (h.row()*P_mat*P_mat.T*a.column())[0,0]/(1+(h.row()*P_mat).norm(p=2)**2)
    return alpha_opt
    
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

# P is a LxL matrix P_mat
def rate_computation_MMSE_alpha(L, M, P_t, H, A):
    P1, P2 = P_t
    # print 'P1 = ', P1, '   P2 = ', P2
    if math.isnan(P1) or math.isnan(P2):
        print 'P1 or P2 should not be NaN!'
        return 0
    P = matrix.diagonal([sqrt(x) for x in P_t])
    
    r = zero_vector(RR, L)
    phi = [float('inf')]*M
    try:
        for i_m in range(0, M):
            h_m = H.row(i_m)
            a_m = A.row(i_m)
            phi[i_m] = (a_m.row()*P).norm(p=2)**2- \
                ((h_m.row()*P*P.T*a_m.column()).norm(p=2)**2)/(1+(h_m*P).norm(p=2)**2)
        for i_l in range(0, L):
            phi_max = 0
            for i_m in range(0, M):
                if A[i_m, i_l] != 0:
                    phi_max = max(phi[i_m], phi_max)
            r[i_l] = 0.5*log(max(1, P[i_l,i_l]**2/phi_max), 2)
    except:
        print 'error in rate_computation_MMSE_alpha'
        raise

    return (r, phi)
    # phi is exactly the fine lattices at the relays

def sum_rate_computation_MMSE_alpha(L, M, P_t, H, A):
    r, phi = rate_computation_MMSE_alpha(L, M, P_t, H, A)
    return sum(r)


# Simple structure, for containing simulation result
class CoF_Sim_Result:
    def __init__(self, sum_rate, sum_rate_var):
        self.sum_rate = sum_rate
        self.sum_rate_var = sum_rate_var
        
class CoF_Dual_Hops_Sim_Result:
    def __init__(self, sum_rate_fixed_pow_sim_mod, sum_rate_sim_mod, sum_rate_opt_mod):
        self.sum_rate_fixed_pow_sim_mod = sum_rate_fixed_pow_sim_mod
        self.sum_rate_sim_mod = sum_rate_sim_mod
        self.sum_rate_opt_mod = sum_rate_opt_mod

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