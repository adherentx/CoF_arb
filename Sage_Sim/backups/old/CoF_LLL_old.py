def CoF_compute_fixed_pow(P_t, is_return_A, *params):
    if len(params) == 2:
        H_a, is_dual_hop = params
    elif len(params) == 4:
        H_a, is_dual_hop, rate_sec_hop, mod_scheme = params
    else:
        raise Exception('error: please check your parameters!')
    
    (M, L) = (H_a.nrows(), H_a.ncols())
    
    try:
        P_t[0]
    except:
        P_t = [P_t]
        
    for i_P in range(0, L):
        if math.isnan(P_t[i_P]):
            print 'P', str(i_P), ' should not be NaN!'
            return 0
        if P_t[i_P] <= 0:
            print 'P', str(i_P), ' should be positive'
            return 0
    
    P_vec = vector(RR, P_t)
    P_mat = matrix.diagonal([sqrt(x) for x in P_vec])
    # Use LLL to find a good A matrix
    # determine the fine lattice of m-th relay at the same time
    try:
        (A_best_LLL, sum_rate_A_LLL, relay_fine_lattices) = Find_A_and_Rate(P_mat, P_vec, H_a, is_return_rate_list=False)
    except:
        print 'error in seeking A and rate'
        raise
    try:
        if is_dual_hop == True:
            '''constraints of the second hop'''
            # relay_fine_lattices is already obtained

            # determine the fine lattice of the l-th transmitter
            trans_fine_lattices = [float(0)]*L
            for i_L in range(0, L):
                for i_M in range(0, M):
                    if (A_best_LLL[i_M, i_L]!=0) and (relay_fine_lattices[i_M]>trans_fine_lattices[i_L]):
                        trans_fine_lattices[i_L] = relay_fine_lattices[i_M]
            # compute the coarse lattice of the l-th transmitter
            trans_coarse_lattices = list(P_vec) # copy
            
            # check whether the second-hop constraint rate_sec_hop can support the first-hop rate r
            try:
                is_support = second_hop(trans_fine_lattices, trans_coarse_lattices, A_best_LLL, rate_sec_hop, mod_scheme)
            except:
                print 'error in second hop'
                raise
        else:
            # if no second hop constraint, then consider the rate r to be supportable
            is_support = True
    except:
        print 'error in checking second hop'
        raise
    
    if is_support == True:
        if is_return_A == True:
            ret = (sum_rate_A_LLL, A_best_LLL)
            return ret
        else:
            return sum_rate_A_LLL

    else:
        if is_return_A == True:
            ret = (0, A_best_LLL)
            return ret
        else:
            return 0

