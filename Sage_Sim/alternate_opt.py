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

def alternate_optimize(P_con, H_a, is_dual_hop, rate_sec_hop=[], mod_scheme='sim_mod'):
    max_iter_alt = 10
    P_vec_init = vector(RR, [P_con, P_con])
    P_mat_init = matrix.diagonal([sqrt(x) for x in P_vec_init])
    try:
        A, sum_rate_0, relay_fine_lattices = Find_A_and_Rate(P_mat_init, P_vec_init, H_a)
        for i_iter in range(0, max_iter_alt):
            P_lst, sum_rate = search_P_for_rate_compute_MMSE_alpha(P_con, A, H_a)
            if sum_rate_0 != sum_rate:
                print 'here'
            P_vec = vector(RR, P_lst)
            P_mat = matrix.diagonal([sqrt(x) for x in P_vec])
            A_old = A
            A, sum_rate, relay_fine_lattices = Find_A_and_Rate(P_mat, P_vec, H_a)
            if A_old == A:
                # converge
                break
    except:
        print 'error in alternate optimize'
        raise
    P_opt = P_vec
    sum_rate_opt = sum_rate
    return sum_rate_opt

def search_P_for_rate_compute_MMSE_alpha(P_con, A, H):
    (M, L) = (H.nrows(), H.ncols())
    cof_pow = lambda x: -sum_rate_computation_MMSE_alpha(L, M, x, H, A)
    Pranges = ((0.1, P_con), (0.1, P_con))
    initial_guess = [0.5*P_con, 0.5*P_con]
    if P_Search_Alg == 'brute':
        res_cof = optimize.brute(cof_pow, Pranges, Ns=20, full_output=True, finish=None)
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










