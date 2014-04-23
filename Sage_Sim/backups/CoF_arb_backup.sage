'''
Simulation of Compute and Forward with arbitrary coarse and fine
lattices.
Author: Yihua Tan
Email: adherentx.tyh@gmail.com
The Chinese University of Hong Kong
'''
import copy
import sys
from sage.all import *
from numpy import arange
import time

print 'Hello, this is the simulation of CoF.'

L = 2 # L transmitters
M = 2 # M relays
p = 7 # The prime number

Cores = 4 # The number of CPU cores used in parallel computing
DEBUG_H = False # When this value is True, the channel matrix H_a is set as certain matrices

def alpha_find(h_m, P_mat, a_m):
    alpha_opt_v = h_m.row()*P_mat*P_mat.T*a_m.column()*~(1+h_m.row()*P_mat*P_mat.T*h_m.column())
    return alpha_opt_v[0, 0]
    
def rate_computation(L, M, P_vec, alpha, H, A):
    r = zero_vector(RR, L)
    for i_l in range(0, L):
        if A.column(i_l).is_zero():
            r[i_l] = 0 # all coefficients are 0.
            raise Exception('It should be impossible')
        else:
            phi_max = 0
            for i_m in range(0, M):
                if A[i_m, i_l] != 0:
                    sum_mis = 0
                    for i_mis in range(0, L):
                        sum_mis = sum_mis+(alpha[i_m]*H[i_m, i_mis]-A[i_m, i_mis])**2*P_vec[i_mis]
                    phi = (alpha[i_m])**2+sum_mis
                    phi_max = max(phi, phi_max)
            r[i_l] = 0.5*log(max(1, P_vec[i_l]/phi_max), 2)
    return r

# Simple structure, for containing simulation result
class CoF_Sim_Result:
    def __init__(self, sum_rate, sum_rate_var):
        self.sum_rate = sum_rate
        self.sum_rate_var = sum_rate_var

class CoF_Sim_Result_for_Fixed_H_Fixed_P:
    def __init__(self, H_a, sum_rate_i_H, A_best, alpha_opt_for_A_best):
        self.H_a = H_a
        self.sum_rate_i_H = sum_rate_i_H
        self.A_best = A_best
        self.alpha_opt_for_A_best = alpha_opt_for_A_best
    def __repr__(self):
        return '---------------Fixed_H_Fixed_P---------------\n' \
            +'H_a = \n'+self.H_a.__str__()+'\n' \
            +'sum_rate_i_H = '+self.sum_rate_i_H.__str__()+'\n' \
            +'A_best = \n'+self.A_best.__str__()+'\n' \
            +'alpha_opt_for_A_best = '+self.alpha_opt_for_A_best.__str__()+'\n' \
            +'--------------------------------------------\n'
        
class CoF_Sim_Result_for_Fixed_H_Variable_P:
    def __init__(self, H_a, sum_rate_i_H_var, A_best_var, \
                 alpha_opt_for_A_best_var, P_vec_best_for_A_best_var):
        self.H_a = H_a
        self.sum_rate_i_H_var = sum_rate_i_H_var
        self.A_best_var = A_best_var
        self.alpha_opt_for_A_best_var = alpha_opt_for_A_best_var
        self.P_vec_best_for_A_best_var = P_vec_best_for_A_best_var

    def __repr__(self):
        return '-------------Fixed_H_Variable_P--------------\n' \
            +'H_a = \n'+self.H_a.__str__()+'\n' \
            +'sum_rate_i_H_var = '+self.sum_rate_i_H_var.__str__()+'\n' \
            +'A_best_var = \n'+ self.A_best_var.__str__()+'\n' \
            +'alpha_opt_for_A_best_var = '+self.alpha_opt_for_A_best_var.__str__()+'\n' \
            +'P_vec_best_for_A_best_var = '+self.P_vec_best_for_A_best_var.__str__()+'\n' \
            +'--------------------------------------------\n'

@parallel(ncpus=Cores)
def CoF_compute_eq_pow_con(P_con):
    iter_H = 20
    division_P = 2 # divided into ? level.
    sum_rate = 0
    sum_rate_var = 0
    for i_H in range(0, iter_H):
        set_random_seed() # to avoid producing the same H_a in different threads
        if DEBUG_H == True:
            H_a = matrix(RR, M, L, [[-0.333414283246484, 0.675839593022605], [0.000374794674703693, 0.766514412738790]])
        else:
            H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
        result_i_H = []
        result_i_H_var = []
        sum_rate_i_H = 0
        sum_rate_i_H_var = 0
        # A_best = zero_matrix(ZZ, M, L)
        # A_best_var = zero_matrix(ZZ, M, L)
        #A_best = []
        #A_best_var = []
        # P_vec_best_for_A_best_var = zero_vector(RR, L)
        #P_vec_best_for_A_best_var = []
        max_a = 2;
        # alpha_opt_for_A_best = zero_vector(RR, M)
        # alpha_opt_for_A_best_var = zero_vector(RR, M)
        #alpha_opt_for_A_best = []
        #alpha_opt_for_A_best_var = []
        for i_A in range(0, (2*max_a+1)**(M*L)):
            k = i_A
            A_temp = zero_matrix(ZZ, M, L)
            for i_M_L in range(0, M*L):
                A_temp[int(i_M_L/M), mod(i_M_L, M)] = int(k/(2*max_a+1)**(M*L-i_M_L-1));
                k = k-A_temp[int(i_M_L/M), mod(i_M_L, M)]*(2*max_a+1)**(M*L-i_M_L-1);
            A = A_temp-max_a*ones_matrix(M, L)
            
            if DEBUG_H == True:
                if i_A == 0:
                    A = matrix(RR, M, L, [[-1, 2], [0, -2]])
                elif i_A == 1:
                    A = matrix(RR, M, L, [[-1, 2], [0, 1]])
            
            A_F = matrix(GF(p), A)
            if A_F.rank() == min(L, M):
                # Fixed power
                P_mat = sqrt(P_con)*identity_matrix(L)
                P_vec = P_con*vector([1]*L)
                
                alpha_opt = zero_vector(RR, M)
                for i_alpha in range(0, M):
                    alpha_opt[i_alpha] = alpha_find(H_a.row(i_alpha), P_mat, A.row(i_alpha))
                r = rate_computation(L, M, P_vec, alpha_opt, H_a, A)
                sum_rate_i_A = sum(r)
                
                P_vec = zero_vector(RR, L)
                P_mat = zero_matrix(RR, L)
                
                # Variable power
                sum_rate_i_A_var = 0;
                P_vec_best = zero_vector(RR, len(P_vec))
                alpha_opt_for_P_best = zero_vector(RR, M)
                delta_P = P_con/(division_P-1)

                
                for i_P_prod in range(0, division_P**L):
                    P_prod = i_P_prod;
                    for i_dim_P in range(0, L):
                        P_prod_t_dim = division_P**(L-i_dim_P-1)
                        P_temp = int(P_prod/P_prod_t_dim)*delta_P
                        P_mat[L-i_dim_P-1, L-i_dim_P-1] = sqrt(P_temp)
                        P_vec[L-i_dim_P-1] = deepcopy(P_temp)
                        P_prod = P_prod - int(P_prod/P_prod_t_dim)*P_prod_t_dim
                    # for i_dim_P
                    alpha_opt_var = zero_vector(RR, M)
                    for i_alpha in range(0, M):
                        alpha_opt_var[i_alpha] = alpha_find(H_a.row(i_alpha), P_mat, A.row(i_alpha))
                    r = rate_computation(L, M, P_vec, alpha_opt_var, H_a, A)
                    if sum_rate_i_A_var < sum(r):
                        sum_rate_i_A_var = sum(r)
                        P_vec_best = deepcopy(P_vec)
                        alpha_opt_for_P_best = deepcopy(alpha_opt_var)
                # for i_P_prod
            # if A_F.is_invertible
            
            else:
                # rank failure: outage
                sum_rate_i_A = 0
                sum_rate_i_A_var = 0
            
            if sum_rate_i_H < sum_rate_i_A:
                result_i_H = [CoF_Sim_Result_for_Fixed_H_Fixed_P(H_a, \
                    sum_rate_i_A, A, alpha_opt)]
                sum_rate_i_H = sum_rate_i_A
                #A_best = [A]
                #alpha_opt_for_A_best = [alpha_opt]
            elif sum_rate_i_H == sum_rate_i_A > 0:
                # Equivalent matrix A
                result_i_H += [CoF_Sim_Result_for_Fixed_H_Fixed_P(H_a, \
                    sum_rate_i_A, A, alpha_opt)]
                sum_rate_i_H = sum_rate_i_A
                #A_best += [A]
                #alpha_opt_for_A_best += [alpha_opt]
            
            if sum_rate_i_H_var < sum_rate_i_A_var:
                result_i_H_var = [CoF_Sim_Result_for_Fixed_H_Variable_P(H_a,\
                    sum_rate_i_A_var, A, alpha_opt_for_P_best, P_vec_best)]
                sum_rate_i_H_var = sum_rate_i_A_var
                #A_best_var = [A]
                #alpha_opt_for_A_best_var = [alpha_opt_for_P_best]
                #P_vec_best_for_A_best_var = [P_vec_best]
            elif sum_rate_i_H_var == sum_rate_i_A_var > 0:
                # Equivalent matrix A
                result_i_H_var += [CoF_Sim_Result_for_Fixed_H_Variable_P(H_a,\
                    sum_rate_i_A_var, A, alpha_opt_for_P_best, P_vec_best)]
                sum_rate_i_H_var = sum_rate_i_A_var
                #A_best_var += [A]
                #alpha_opt_for_A_best_var += [alpha_opt_for_P_best]
                #P_vec_best_for_A_best_var += [P_vec_best]
            
        # i_A
        '''
        print '*******************************************************'
        print 'Power Constraint = '; print P_con*vector([1]*L)
        print 'H_a = '; print H_a
        print '-------------------------------------------------------'
        print 'alpha_opt_for_A_best = '; print alpha_opt_for_A_best
        print 'A_best = '; print A_best
        print 'sum_rate_i_H = '; print sum_rate_i_H
        print '-------------------------------------------------------'
        print 'alpha_opt_for_A_best_var = '; print alpha_opt_for_A_best_var
        print 'A_best_var = '; print A_best_var
        print 'P_vec_best_for_A_best_var = '; print P_vec_best_for_A_best_var
        print 'sum_rate_i_H_var = '; print sum_rate_i_H_var
        '''
        print result_i_H
        print result_i_H_var
        
        sum_rate += sum_rate_i_H
        sum_rate_var += sum_rate_i_H_var
    # for i_H
    sum_rate /= iter_H
    sum_rate_var /= iter_H
    
    return CoF_Sim_Result(sum_rate, sum_rate_var)


P_eq_dB_Min = float(20)
P_eq_dB_Max = float(40)
P_delta = 5
P_eq_dB = arange(P_eq_dB_Min, P_eq_dB_Max, P_delta)
P_eq = [10**(P_eq_dB_i/10) for P_eq_dB_i in P_eq_dB]
Pl_con = P_eq

t1 = time.ctime()
result = list(CoF_compute_eq_pow_con(Pl_con))
t2 = time.ctime()

print 'Simulation started at ', t1
print 'Simulation ended at  ', t2

sum_rate = [result[i][1].sum_rate for i in range(0, len(Pl_con))]
sum_rate_var = [result[i][1].sum_rate_var for i in range(0, len(Pl_con))]
        
plot_sum_rate = list_plot(zip(P_eq_dB, sum_rate), plotjoined=True, marker='o', \
                          rgbcolor=Color('red'), linestyle="--")
plot_sum_rate_var = list_plot(zip(P_eq_dB, sum_rate_var), plotjoined=True, marker='x', \
                              rgbcolor=Color('blue'), linestyle='-')
plot_compare = plot_sum_rate+plot_sum_rate_var
str_label = t1
plot_compare.save('/home/adherentx/Dropbox/Research/My_Report/Compute_and_Forward/Sage_Sim/Comparison_Fixed_and_Variable_Power' \
                  +str_label+'.eps')
show(plot_compare)



