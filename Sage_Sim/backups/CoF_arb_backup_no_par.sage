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

print 'Hello, this is the simulation of CoF.'

L = 2 # L transmitters
M = 2 # M relays
p = 7 # The prime number

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
                    phi = (10*alpha[i_m])**2+sum_mis
                    phi_max = max(phi, phi_max)
            r[i_l] = 0.5*log(max(1, P_vec[i_l]/phi_max), 2)
    return r
    
P_eq_dB_Min = 20
P_eq_dB_Max = 30
P_delta = 5
P_eq_dB = range(P_eq_dB_Min, P_eq_dB_Max, P_delta)
P_eq = [10**(P_eq_dB_i/10) for P_eq_dB_i in P_eq_dB]
Pl_con = P_eq

iter_H = 500
sum_rate_P = zero_vector(RR, 3)

sum_rate = zero_vector(RR, len(P_eq))
sum_rate_var = zero_vector(RR, len(P_eq))

for i_P in range(0, len(P_eq)):
    print 'Simulating ', i_P, '/', len(P_eq)
    for i_H in range(0, iter_H):
        # H_a = matrix([[1.5, -1], [2, 1]])
        H_a = matrix.random(RR, M, L)
        sum_rate_i_H = 0
        sum_rate_i_H_var = 0
        A_best = zero_matrix(ZZ, M, L)
        A_best_var = zero_matrix(ZZ, M, L)
        P_vec_best_for_A_best_var = zero_vector(RR, L)
        max_a = 4;
        alpha_opt_for_A_best = zero_vector(RR, M)
        alpha_opt_for_A_best_var = zero_vector(RR, M)
        for i_A in range(0, (2*max_a+1)**(M*L)):
            k = i_A
            A_temp = zero_matrix(ZZ, M, L)
            for i_M_L in range(0, M*L):
                A_temp[int(i_M_L/M), mod(i_M_L, M)] = int(k/(2*max_a+1)**(M*L-i_M_L-1));
                k = k-A_temp[int(i_M_L/M), mod(i_M_L, M)]*(2*max_a+1)**(M*L-i_M_L-1);
            A = A_temp-max_a*ones_matrix(M, L)
            
            if A == matrix([[-1, 0], [-2, -1]]):
                pass
            if A == matrix([[3, -2], [2, 1]]):
                pass
            
            A_F = matrix(GF(p), A)
            if A_F.is_invertible():
                P_vec = zero_vector(RR, L)
                P_mat = zero_matrix(RR, L)
                # Fixed power
                
                P_mat = sqrt(Pl_con[i_P])*identity_matrix(L)
                P_vec = Pl_con[i_P]*vector([1]*L)
                
                alpha_opt = zero_vector(RR, M)
                for i_alpha in range(0, M):
                    alpha_opt[i_alpha] = alpha_find(H_a.row(i_alpha), P_mat, A.row(i_alpha))
                #r = zero_vector(RR, L)
                r = rate_computation(L, M, P_vec, alpha_opt, H_a, A)
                sum_rate_i_A = sum(r)
                
                
                P_vec = zero_vector(RR, L)
                P_mat = zero_matrix(RR, L)
                
                # Variable power
                
                sum_rate_i_A_var = 0;
                P_vec_best = zero_vector(RR, len(P_vec))
                alpha_opt_for_P_best = zero_vector(RR, M)
                division_P = 2 # divided into ? level.
                delta_P = Pl_con[i_P]/(division_P-1)
                for i_P_prod in range(0, division_P**L):
                    P_prod = i_P_prod;
                    for i_dim_P in range(0, L):
                        P_prod_t_dim = division_P**(L-i_dim_P-1)
                        P_temp = int(P_prod/P_prod_t_dim)*delta_P
                        P_mat[L-i_dim_P-1, L-i_dim_P-1] = sqrt(P_temp)
                        P_vec[L-i_dim_P-1] = deepcopy(P_temp)
                        P_prod = P_prod - int(P_prod/P_prod_t_dim)*P_prod_t_dim
                    # for i_dim_P
                    # print P_vec
                    if P_vec == vector([10000, 10000]) and A == matrix([[1, 0], [2, 1]]):
                        pass
                    if A == matrix([[0, -4], [-1, -4]]):
                        pass
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
                sum_rate_i_H = sum_rate_i_A
                A_best = A
                alpha_opt_for_A_best = alpha_opt
            
            
            if sum_rate_i_H_var < sum_rate_i_A_var:
                sum_rate_i_H_var = sum_rate_i_A_var
                A_best_var = A
                alpha_opt_for_A_best_var = alpha_opt_for_P_best
                P_vec_best_for_A_best_var = P_vec_best
            
        # i_A
        
        print 'Power Constraint = '; print Pl_con[i_P]*vector([1]*L)
        print 'H_a = '; print H_a
        
        
        print 'alpha_opt_for_A_best = '; print alpha_opt_for_A_best
        print 'A_best = '; print A_best
        
        
        print 'alpha_opt_for_A_best_var = '; print alpha_opt_for_A_best_var
        print 'A_best_var = '; print A_best_var
        print 'P_vec_best_for_A_best_var = '; print P_vec_best_for_A_best_var
        print 'sum_rate_i_H_var = '; print sum_rate_i_H_var
        
        sum_rate[i_P] += sum_rate_i_H
        sum_rate_var[i_P] += sum_rate_i_H_var
    # for i_H
# for i_P
sum_rate /= iter_H
sum_rate_var /= iter_H
        
plot_sum_rate = list_plot(zip(P_eq_dB, sum_rate), plotjoined=True, marker='o', \
                          rgbcolor=Color('red'), linestyle="--")
plot_sum_rate_var = list_plot(zip(P_eq_dB, sum_rate_var), plotjoined=True, marker='x', \
                              rgbcolor=Color('blue'), linestyle='-')
show(plot_sum_rate+plot_sum_rate_var)



