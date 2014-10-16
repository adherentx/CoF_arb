
def CoF_compute_search_pow(P_con, H_a, is_dual_hop, rate_sec_hop=[], mod_scheme='sim_mod'):
    (M, L) = (H_a.nrows(), H_a.ncols())
    cof_pow = lambda x: -CoF_compute_fixed_pow(x, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme)
    Pranges = ((0.1, P_con), )*L
    initial_guess = [0.5*P_con]*L
    try:
        if P_Search_Alg == 'brute':
            res_cof = optimize.brute(cof_pow, Pranges, Ns=brute_number, full_output=True, finish=None)
            P_opt = res_cof[0]
            sum_rate_opt = -res_cof[1] # negative! see minus sign in cof_pow
        elif P_Search_Alg == 'TNC':
            #res_cof = optimize.minimize(cof_pow, initial_guess, method='TNC', bounds=Pranges, options={'maxiter': 400, 'approx_grad': True})
            #P_opt = list(res_cof.x)
            #sum_rate_opt = -res_cof.fun # negative! see minus sign in cof_pow
            res_cof = optimize.fmin_tnc(cof_pow, initial_guess, bounds=list(Pranges), approx_grad=True, epsilon=1, stepmx=10)
            P_opt = res_cof[0]
            sum_rate_opt = CoF_compute_fixed_pow(P_opt, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme)
        elif P_Search_Alg == 'anneal':
            res_cof = optimize.anneal(cof_pow, initial_guess, schedule='cauchy', T0=1, Tf=1e-6, \
                      full_output=True, maxiter=30, lower=[1, 1], upper=[P_con, P_con], dwell=30, disp=True)
            P_opt = list(res_cof[0])
            sum_rate_opt = -res_cof[1]
        else:
            raise Exception('error: algorithm not supported')
    except:
        print 'error in search algorithms'
        raise
    return sum_rate_opt