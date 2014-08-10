'''
This includes some functions related with the second hop in CoF.
'''
from sage.all import *
import copy
from numpy import arange
from sage.parallel.all import *
import time
from CoF_basic import *
import itertools



def second_hop_support_rates(relay_fine_lattices, trans_coarse_lattices, A, rate_sec_hop, mod_scheme):
    (M, L) = (A.nrows(), A.ncols())
    if M != L:
        raise Exception("L and M should be the same in destination's perspective.")
    
    if mod_scheme == 'opt_mod':
        # iterate all permutations A^L_L. 
        sum_rate_max = 0
        has_feasible = False
        for mod_order in itertools.permutations(list(range(0, M)), L):
            relay_actual_fine_lattices = list(relay_fine_lattices)
            # each element in mod_order: the relay that the l-th lattice should be assigned to
            mod_order_list = list(mod_order) # note: not equal to pi_e()
            is_decodable = True
            A_F = matrix(GF(p), A)
            A_i_L = A_F
            x = list(trans_coarse_lattices) # copy
            for i_L in range(0, L):
                if A_i_L.rank() < L-i_L:
                    is_decodable = False
                    break
                elif A_i_L.rank() == L-i_L:
                    x_min_idx, x_min_val = min(enumerate(x), key=lambda x:x[1])
                    x.pop(x_min_idx)
                    if A_i_L.rank() > 1:
                        A_i_L = A_i_L.delete_columns([x_min_idx])
                        row_thresh = mod_order_list[x_min_idx]
                        A_i_L = A_i_L.delete_rows([mod_order_list.pop(x_min_idx)])
                        for i_row in range(0, len(mod_order_list)):
                            if mod_order_list[i_row] > row_thresh:
                                mod_order_list[i_row] -= 1
                            pass
                    else:
                        # don't need to reduce since no further operations need it
                        pass
                else:
                    raise Exception('error: the rank should not exceed L-i_L')
            if is_decodable == True:
                has_feasible = True
                # if decodable, then calculate how much sum rate it can support
                mod_order_list = list(mod_order)
                # determine the coarse lattices at the relays
                relay_coarse_lattices = [0]*M
                for i_m in range(0, M):
                    relay_coarse_lattices[i_m] = trans_coarse_lattices[mod_order_list.index(i_m)]
                # determine the achievable fine lattices at the relays
                for i_m in range(0, M):
                    R_m = max(0, 0.5*log(relay_coarse_lattices[i_m]/relay_actual_fine_lattices[i_m], 2))
                    if rate_sec_hop[i_m] < R_m:
                        relay_actual_fine_lattices[i_m] = relay_coarse_lattices[i_m]/(2**(2*rate_sec_hop[i_m]))
                # determine the fine lattice of the l-th transmitter
                trans_fine_lattices = [float(0)]*L
                for i_L in range(0, L):
                    for i_M in range(0, M):
                        if (A[i_M, i_L]!=0) and (relay_actual_fine_lattices[i_M]>trans_fine_lattices[i_L]):
                            trans_fine_lattices[i_L] = relay_actual_fine_lattices[i_M]
                # calculate transmission sum rates
                r = [0]*L
                for i_l in range(0, L):
                    r[i_l] = max(0, 0.5*log(trans_coarse_lattices[i_l]/trans_fine_lattices[i_l], 2))
                sum_rate = sum(r)
                if sum_rate_max < sum_rate:
                    sum_rate_max = sum_rate
        # if no successful scheme found,
        if has_feasible == False: 
            error('If Q is full rank, then there must be at leat one feasible scheme. But no one found.')
        return sum_rate_max
    elif mod_scheme == 'naive_mod':
        relay_actual_fine_lattices = list(relay_fine_lattices)
        relay_coarse_lattice = max(trans_coarse_lattices)
        # determine the achievable fine lattices at the relays
        for i_m in range(0, M):
            R_m = max(0, 0.5*log(relay_coarse_lattice/relay_actual_fine_lattices[i_m], 2))
            if rate_sec_hop[i_m] < R_m:
                relay_actual_fine_lattices[i_m] = relay_coarse_lattice/(2**(2*rate_sec_hop[i_m]))
        # determine the fine lattice of the l-th transmitter
        trans_fine_lattices = [float(0)]*L
        for i_L in range(0, L):
            for i_M in range(0, M):
                if (A[i_M, i_L]!=0) and (relay_actual_fine_lattices[i_M]>trans_fine_lattices[i_L]):
                    trans_fine_lattices[i_L] = relay_actual_fine_lattices[i_M]
        # calculate transmission sum rates
        r = [0]*L
        for i_l in range(0, L):
            r[i_l] = max(0, 0.5*log(trans_coarse_lattices[i_l]/trans_fine_lattices[i_l], 2))
        sum_rate = sum(r)
        return sum_rate
    else:
        raise Exception("mod_scheme should take value as 'opt_mod' or 'sim_opt'!")
    

def sum_rate_computation_MMSE_alpha_two_hop(L, M, P_t, H, A, is_dual_hop, rate_sec_hop, mod_scheme):
    r, relay_fine_lattices = rate_computation_MMSE_alpha(L, M, P_t, H, A)
    
    '''constraints of the second hop'''
    # relay_fine_lattices is already obtained

    # determine the fine lattice of the l-th transmitter
    trans_fine_lattices = [float(0)]*L
    for i_L in range(0, L):
        for i_M in range(0, M):
            if (A[i_M, i_L]!=0) and (relay_fine_lattices[i_M]>trans_fine_lattices[i_L]):
                trans_fine_lattices[i_L] = relay_fine_lattices[i_M]
    # compute the coarse lattice of the l-th transmitter
    trans_coarse_lattices = list(P_t) # copy
    
    if is_dual_hop == True:
        # check whether the second-hop constraint rate_sec_hop can support the first-hop rate r
        try:
            is_support = second_hop(trans_fine_lattices, trans_coarse_lattices, A, rate_sec_hop, mod_scheme)
        except:
            print 'error in second hop. in sum_rate_computation_MMSE_alpha_two_hop'
            raise
    else:
        is_support = True
    
    if is_support == True:
        return sum(r)
    else:
        return 0

            

if __name__ == "__main__":
    print '-----------------------------------\n'+ \
        'testing CoF_second_hop\n'
#     fine_lattices = [0.5, 0.2]
#     coarse_lattices = [1.5, 2.1]
#     R = [1.5, 1.7]
#     A = matrix(ZZ, 2, 2, [[1, 2], [1, 1]])
#     p = 3 
#     print second_hop(fine_lattices, coarse_lattices, A, R, 'opt_mod')
#     
#     fine_lattices = [0.5, 0.2, 0.4]
#     coarse_lattices = [1.5, 2.1, 1]
#     R = [1.5, 1.7, 3]
#     A = matrix(ZZ, 3, 3, [[0, 0, 2], [2, 2, 0], [0, 2, 2]])
#     print second_hop(fine_lattices, coarse_lattices, A, R, 'opt_mod')
    
    R = [1, 2]
    A = matrix(ZZ, 2, 2, [[1, 2], [1, 1]])
    p = 3
    relay_fine_lattices = [0.5, 0.2]
    trans_coarse_lattices = [1.5, 2.0]
    print RR(second_hop_support_rates(relay_fine_lattices, trans_coarse_lattices, A, R, 'opt_mod'))
    # should be 1.792
    
    R = [1, 10]
    A = matrix(ZZ, 2, 2, [[1, 2], [0, 1]])
    p = 3
    relay_fine_lattices = [0.5, 0.2]
    trans_coarse_lattices = [1.5, 2.0]
    print RR(second_hop_support_rates(relay_fine_lattices, trans_coarse_lattices, A, R, 'opt_mod'))
    # should be 1.792

    R = [1, 1]
    A = matrix(ZZ, 2, 2, [[0, 2], [1, 1]])
    p = 3
    relay_fine_lattices = [0.5, 0.2]
    trans_coarse_lattices = [1.5, 2.0]
    print RR(second_hop_support_rates(relay_fine_lattices, trans_coarse_lattices, A, R, 'opt_mod'))
    # should be 1.792
    