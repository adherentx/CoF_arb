'''
This is something related with alternating optimization.
'''
from sage.all import *
import copy
from numpy import arange
from sage.parallel.all import *
import time
from CoF_basic import *
from CoF_LLL import *
from scipy import optimize


def alternate_optimize(P_con, H_a, is_dual_hop, rate_sec_hop=[], mod_scheme='sim_mod', quan_scheme='sym_quan'):
    (M, L) = (H_a.nrows(), H_a.ncols())
    max_iter_alt = 10
    # P_t_init = [0.8*P_con]*L
    h_pow = [norm(H_a.column(i_col))**2 for i_col in range(0, L)]
    h_pow_max = max(h_pow)
    P_t_init = [h_pow[i]/h_pow_max*P_con for i in range(0, L)]
    sum_rate = 0
    try:
        P_lst = P_t_init
        A_old = 0 # deliberately. to avoid the problem of A=zero_matrix at the first iteration
        for i_iter in range(0, max_iter_alt):
            sum_rate_A, A = CoF_compute_fixed_pow(P_lst, True, H_a, is_dual_hop, rate_sec_hop, mod_scheme)
            P_lst, sum_rate_P = search_P_for_rate_compute_MMSE_alpha(P_con, A, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme)
            if A_old == A:
                # converge
                break
            A_old = A
            if i_iter == max_iter_alt:
                raise("It doesn't converge in "+str(max_iter_alt)+' iterations!')
    except:
        print 'error in alternate optimize'
        raise
    P_opt = P_lst
    sum_rate_opt = sum_rate_P
    return sum_rate_opt

# FIX Me!!!
def search_P_for_rate_compute_MMSE_alpha(P_con, A, H, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme='sym_quan'):
    (M, L) = (H.nrows(), H.ncols())
    cof_pow = lambda x: -sum_rate_computation_MMSE_alpha_two_hop(L, M, x, H, A, is_dual_hop, rate_sec_hop, mod_scheme)
    Pranges = ((0.1, P_con), )*L
    initial_guess = [0.5*P_con]*L
    if P_Search_Alg == 'brute':
        res_cof = optimize.brute(cof_pow, Pranges, Ns=brute_number, full_output=True, finish=None)
        P_opt = list(res_cof[0])
        sum_rate_opt = -res_cof[1] # negative! see minus sign in cof_pow
    elif P_Search_Alg == 'TNC':
        res_cof = optimize.minimize(cof_pow, initial_guess, method='TNC', bounds=Pranges, options={'maxiter': 100})
        P_opt = list(res_cof.x)
        sum_rate_opt = -res_cof.fun # negative! see minus sign in cof_pow
    elif P_Search_Alg == 'anneal':
        res_cof = optimize.anneal(cof_pow, initial_guess, schedule='boltzmann', \
                                  full_output=True, maxiter=20, lower=2, upper=P_con, dwell=20, disp=True)
        P_opt = list(res_cof[0])
        sum_rate_opt = -res_cof[1]
    else:
        raise Exception('error: algorithm not supported')
    return (P_opt, sum_rate_opt)










