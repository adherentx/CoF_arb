'''
Simulation of Compute and Forward with arbitrary coarse and fine
lattices.
Author: Yihua Tan
Email: adherentx.tyh@gmail.com
The Chinese University of Hong Kong
'''


''' MEMO:
  Fix all the 'FIX ME' bugs!!!
'''


import copy
import sys
from sage.all import *
from numpy import arange
import time
from CoF_LLL import *
from CoF_basic import *
from CoF_second_hop import *
from scipy import optimize
import math

print 'Hello, this is the simulation of CoF.'

def CoF_compute_search_pow(P_con, H_a, is_dual_hop, rate_sec_hop=[], mod_scheme='sim_mod'):
    cof_pow = lambda x: -CoF_compute_fixed_pow(x, H_a, is_dual_hop, rate_sec_hop, mod_scheme)
    Pranges = ((0.1, P_con), (0.1, P_con))
    initial_guess = [0.5*P_con, 0.5*P_con]
    if P_Search_Alg == 'brute':
        res_cof = optimize.brute(cof_pow, Pranges, Ns=20, full_output=True, finish=None)
        P_opt = res_cof[0]
        sum_rate_opt = -res_cof[1] # negative! see minus sign in cof_pow
    elif P_Search_Alg == 'TNC':
        res_cof = optimize.minimize(cof_pow, initial_guess, method='TNC', bounds=Pranges, options={'maxiter': 100})
        P_opt = list(res_cof.x)
        sum_rate_opt = -res_cof.fun # negative! see minus sign in cof_pow
    elif P_Search_Alg == 'anneal':
        res_cof = optimize.anneal(cof_pow, initial_guess, schedule='boltzmann', full_output=True, maxiter=20, lower=2, upper=P_con, dwell=20, disp=True)
        P_opt = list(res_cof[0])
        sum_rate_opt = -res_cof[1]
    else:
        raise Exception('error: algorithm not supported')
    return sum_rate_opt
    
    
'''This is for L=M=2!'''
def CoF_compute_fixed_pow(P_t, *params):
    P1, P2 = P_t
    print 'P1 = ', P1, '   P2 = ', P2
    if math.isnan(P1) or math.isnan(P2):
        print 'P1 or P2 should not be NaN!'
        return 0

    if len(params) == 2:
        H_a, is_dual_hop = params
    elif len(params) == 4:
        H_a, is_dual_hop, rate_sec_hop, mod_scheme = params
    else:
        raise Exception('error: please check your parameters!')
    
    if P1 <= 0 or P2 <= 0:
        print 'P1 and P2 should be positive'
        return 0
        #raise Exception('error: P1 and P2 should be positive')
    
    P_vec = vector(RR, [P1, P2])
    P_mat = matrix.diagonal([sqrt(x) for x in P_vec])
    # Use LLL to find a good A matrix
    (A_best_LLL, sum_rate_A_LLL, alpha_opt_LLL) = Find_A_and_Rate(P_mat, P_vec, H_a)
    if is_dual_hop == True:
        '''constraints of the second hop'''
        # determine the fine lattice of m-th relay
        relay_fine_lattices = [float('inf')]*M
        for i_M in range(0, M):
            relay_fine_lattices[i_M] = alpha_opt_LLL[i_M]**2+\
                (((alpha_opt_LLL[i_M]*H_a.row(i_M)-A_best_LLL.row(i_M))*P_mat).norm())**2
        # determine the fine lattice of the l-th transmitter
        trans_fine_lattices = [float(0)]*L
        for i_L in range(0, L):
            for i_M in range(0, M):
                if (A_best_LLL[i_M, i_L]!=0) and (relay_fine_lattices[i_M]>trans_fine_lattices[i_L]):
                    trans_fine_lattices[i_L] = relay_fine_lattices[i_M]
        # compute the coarse lattice of the l-th transmitter
        trans_coarse_lattices = list(P_vec) # copy
        
        # check whether the second-hop constraint rate_sec_hop can support the first-hop rate r
        is_support = second_hop(trans_fine_lattices, trans_coarse_lattices, A_best_LLL, rate_sec_hop, mod_scheme)
    else:
        # if no second hop constraint, then consider the rate r to be supportable
        is_support = True
    
    if is_support == True:
        # return (A_best_LLL, sum_rate_A_LLL, alpha_opt_LLL)
        return sum_rate_A_LLL
    else:
        # return (zero_matrix(ZZ, M, L), 0, zero_vector(RR, M))
        return 0

@parallel(ncpus=Cores)
def CoF_compute_eq_pow_con_first_hop(P_con):
    iter_H = 2000
    sum_rate = 0
    sum_rate_var = 0
    for i_H in range(0, iter_H):
        set_random_seed() # to avoid producing the same H_a in different threads
        if DEBUG_H == True:
            H_a = matrix(RR, M, L, [[0.430020183081381, 0.507485287402749], [0.571045084182643, -0.846256047736565]])
        else:
            H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
        
        sum_rate_i_H = 0
        sum_rate_i_H_var = 0

        # Fixed power
#         P_mat = sqrt(P_con)*identity_matrix(L)
#         P_vec = P_con*vector([1]*L)
        # (A, sum_rate_A, alpha_opt) = CoF_compute_fixed_pow(P_mat, P_vec, H_a, is_dual_hop=False)
        sum_rate_A = CoF_compute_fixed_pow((P_con, P_con), H_a, False)
        
#         P_vec = zero_vector(RR, L)
#         P_mat = zero_matrix(RR, L)
        
        # Variable power
        sum_rate_A_var = CoF_compute_search_pow(P_con, H_a, is_dual_hop=False)
        '''
        sum_rate_A_var = 0;
        P_vec_best = zero_vector(RR, len(P_vec))
        alpha_opt_for_P_best = zero_vector(RR, M)
        delta_P = P_con/(division_P-1)
        
        for i_P_prod in range(0, division_P**L):
            P_prod = i_P_prod;
            for i_dim_P in range(0, L):
                P_prod_t_dim = division_P**(L-i_dim_P-1)
                P_temp = int(P_prod/P_prod_t_dim)*delta_P
                # FIX ME!
                # delete +0.01 in the following two lines!!!!
                P_mat[L-i_dim_P-1, L-i_dim_P-1] = sqrt(P_temp+0.0001)
                P_vec[L-i_dim_P-1] = deepcopy(P_temp+0.0001)
                P_prod = P_prod - int(P_prod/P_prod_t_dim)*P_prod_t_dim
            # for i_dim_P

            (A_best_LLL, sum_rate_A_LLL, alpha_opt_LLL) = CoF_compute_fixed_pow(P_mat, P_vec, H_a, is_dual_hop=False)
            
            if sum_rate_A_var < sum_rate_A_LLL:
                sum_rate_A_var = sum_rate_A_LLL
                P_vec_best = deepcopy(P_vec)
                alpha_opt_for_P_best = deepcopy(alpha_opt_LLL)
                A_var = A_best_LLL
                    
        # for i_P_prod
        '''
        
        
#         result_i_H = CoF_Sim_Result_for_Fixed_H_Fixed_P(H_a, \
#             sum_rate_A, A, alpha_opt, P_con*vector([1]*L))
        sum_rate_i_H = sum_rate_A
        
#         result_i_H_var = CoF_Sim_Result_for_Fixed_H_Variable_P(H_a,\
#             sum_rate_A_var, A_var, alpha_opt_for_P_best, P_vec_best)
        sum_rate_i_H_var = sum_rate_A_var
        
#         print result_i_H
#         print result_i_H_var
        
        sum_rate += sum_rate_i_H
        sum_rate_var += sum_rate_i_H_var
    # for i_H
    sum_rate /= iter_H
    sum_rate_var /= iter_H
    
    return CoF_Sim_Result(sum_rate, sum_rate_var)

@parallel(ncpus=Cores)
def CoF_compute_eq_pow_con_dual_hops(P_con):
    iter_H = 20000
    sum_rate_sim_mod = 0
    sum_rate_opt_mod = 0
    '''How to determine the power of the relays?'''
    P_relay = 0.25*P_con 
    for i_H in range(0, iter_H):
        set_random_seed() # to avoid producing the same H_a in different threads
        if DEBUG_H == True:
            H_a = matrix(RR, M, L, [[-0.333414283246484, 0.675839593022605], [0.000374794674703693, 0.766514412738790]])
        else:
            H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
        
        '''Fix Me: set the rate constraint of the second hop'''
        set_random_seed()
        H_b = (matrix.random(RR, M, 1, distribution=RealDistribution('gaussian', 1))).column(0)
        rate_sec_hop = [0]*M # [4, 6] # ? bit/s for each parallel channel in the second hop
        for i_h_b in range(0, M):
            rate_sec_hop[i_h_b] = 0.5*log(1+H_b[i_h_b]**2*P_relay, 2)
        
#         P_vec = zero_vector(RR, L)
#         P_mat = zero_matrix(RR, L)
        
        # In dual-hop system, we always use variable power method to avoid outage at the second hop.
        sum_rate_i_H_sim_mod = CoF_compute_search_pow(P_con, H_a, is_dual_hop=True, rate_sec_hop=rate_sec_hop, mod_scheme='sim_mod')
        sum_rate_i_H_opt_mod = CoF_compute_search_pow(P_con, H_a, is_dual_hop=True, rate_sec_hop=rate_sec_hop, mod_scheme='opt_mod')
        '''
        sum_rate_i_H_opt_mod = 0;
        sum_rate_i_H_sim_mod = 0;
        A_opt_mod = zero_matrix(ZZ, 3, 3)
        A_sim_mod = zero_matrix(ZZ, 3, 3)
        P_vec_best_opt_mod = zero_vector(RR, len(P_vec))
        P_vec_best_sim_mod = zero_vector(RR, len(P_vec))
        alpha_opt_opt_mod = zero_vector(RR, M)
        alpha_opt_sim_mod = zero_vector(RR, M)
        delta_P = P_con/(division_P-1)
        
        for i_P_prod in range(0, division_P**L):
            P_prod = i_P_prod;
            for i_dim_P in range(0, L):
                P_prod_t_dim = division_P**(L-i_dim_P-1)
                P_temp = int(P_prod/P_prod_t_dim)*delta_P
                # FIX ME!
                # delete +0.01 in the following two lines!!!!
                P_mat[L-i_dim_P-1, L-i_dim_P-1] = sqrt(P_temp+0.0001)
                P_vec[L-i_dim_P-1] = deepcopy(P_temp+0.0001)
                P_prod = P_prod - int(P_prod/P_prod_t_dim)*P_prod_t_dim
            # for i_dim_P
            
            (A_i_P_opt_mod, sum_rate_i_P_opt_mod, alpha_opt_i_P_opt_mod) = \
                CoF_compute_fixed_pow(P_mat, P_vec, H_a, is_dual_hop=True, \
                rate_sec_hop=rate_sec_hop, mod_scheme='opt_mod')
            (A_i_P_sim_mod, sum_rate_i_P_sim_mod, alpha_opt_i_P_sim_mod) = \
                CoF_compute_fixed_pow(P_mat, P_vec, H_a, is_dual_hop=True, \
                rate_sec_hop=rate_sec_hop, mod_scheme='sim_mod')
            
            if sum_rate_i_H_opt_mod < sum_rate_i_P_opt_mod:
                sum_rate_i_H_opt_mod = sum_rate_i_P_opt_mod
                P_vec_best_opt_mod = deepcopy(P_vec)
                alpha_opt_opt_mod = deepcopy(alpha_opt_i_P_opt_mod)
                A_opt_mod = A_i_P_opt_mod
            if sum_rate_i_H_sim_mod < sum_rate_i_P_sim_mod:
                sum_rate_i_H_sim_mod = sum_rate_i_P_sim_mod
                P_vec_best_sim_mod = deepcopy(P_vec)
                alpha_opt_sim_mod = deepcopy(alpha_opt_i_P_sim_mod)
                A_sim_mod = A_i_P_sim_mod
                    
        # for i_P_prod
        '''

#         result_i_H_opt_mod = CoF_Sim_Result_for_Fixed_H_Variable_P(H_a,\
#             sum_rate_i_H_opt_mod, A_opt_mod, alpha_opt_opt_mod, P_vec_best_opt_mod)
#         result_i_H_sim_mod = CoF_Sim_Result_for_Fixed_H_Variable_P(H_a,\
#             sum_rate_i_H_sim_mod, A_sim_mod, alpha_opt_sim_mod, P_vec_best_sim_mod)
        
#         print result_i_H_opt_mod
#         print result_i_H_sim_mod
        
        sum_rate_sim_mod += sum_rate_i_H_sim_mod
        sum_rate_opt_mod += sum_rate_i_H_opt_mod
    # for i_H
    sum_rate_sim_mod /= iter_H
    sum_rate_opt_mod /= iter_H
    
    return CoF_Sim_Result(sum_rate_sim_mod, sum_rate_opt_mod)


if __name__ == "__main__": 
    '''Equal Power Constraint'''
    P_eq_dB_Min = float(20)
    P_eq_dB_Max = float(60)
    P_delta = 5
    P_eq_dB = arange(P_eq_dB_Min, P_eq_dB_Max, P_delta)
    P_eq = [10**(P_eq_dB_i/10) for P_eq_dB_i in P_eq_dB]
    Pl_con = P_eq
    
    '''First Hop'''
    if True:
        t1 = time.ctime()
        result = list(CoF_compute_eq_pow_con_first_hop(Pl_con))
        t2 = time.ctime()
        
        print 'Simulation of the first hop started at ', t1
        print 'Simulation of the first hop ended at ', t2
        
        sum_rate = [result[i][1].sum_rate_1 for i in range(0, len(Pl_con))]
        sum_rate_var = [result[i][1].sum_rate_2 for i in range(0, len(Pl_con))]
                
        plot_sum_rate = list_plot(zip(P_eq_dB, sum_rate), plotjoined=True, marker='o', \
                                  rgbcolor=Color('red'), linestyle="--", \
                                  legend_label= 'sum rate of fixed power method', \
                                  title = 'Comparison of fixed and variable methods int the first hop')
        plot_sum_rate_var = list_plot(zip(P_eq_dB, sum_rate_var), plotjoined=True, marker='x', \
                                      rgbcolor=Color('blue'), linestyle='-', \
                                      legend_label = 'sum rate of variable power method')
        plot_compare = plot_sum_rate+plot_sum_rate_var
        plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
        str_label = t1
        plot_compare.save('/home/adherentx/Dropbox/Research/My_Report/Compute_and_Forward/Sage_Sim/Simulation_Results/Comparison_Fixed_and_Variable_Power_in_the_First_Hop-' \
                          +P_Search_Alg+str_label+'.eps')
        show(plot_compare)
    
    '''Dual Hops'''
    if False:
        t1 = time.ctime()
        result = list(CoF_compute_eq_pow_con_dual_hops(Pl_con))
        t2 = time.ctime()
        
        print 'Simulation of dual hops started at ', t1
        print 'Simulation of dual hops ended at ', t2
        
        sum_rate_sim_mod = [result[i][1].sum_rate_1 for i in range(0, len(Pl_con))]
        sum_rate_opt_mod = [result[i][1].sum_rate_2 for i in range(0, len(Pl_con))]
                
        plot_sum_rate_sim_mod = list_plot(zip(P_eq_dB, sum_rate_sim_mod), plotjoined=True, marker='o', \
                                  rgbcolor=Color('red'), linestyle="--", \
                                  legend_label= 'sum rate of simple modulo method', \
                                  title = 'Comparison of simple and optimal modulo methods in dual-hops system')
        plot_sum_rate_opt_mod = list_plot(zip(P_eq_dB, sum_rate_opt_mod), plotjoined=True, marker='x', \
                                      rgbcolor=Color('blue'), linestyle='-', \
                                      legend_label = 'sum rate of optimal mudulo method')
        plot_compare = plot_sum_rate_sim_mod+plot_sum_rate_opt_mod
        plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
        str_label = t1
        plot_compare.save('/home/adherentx/Dropbox/Research/My_Report/Compute_and_Forward/Sage_Sim/Simulation_Results/Comparison_Simple_and_Optimal_Modulo_Methods_in_Dual_Hops_System-' \
                          +P_Search_Alg+str_label+'.eps')
        show(plot_compare)
    
    '''Dual Hop'''
    
    raw_input() # stop Sage from shutting down
