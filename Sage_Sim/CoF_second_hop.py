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

# coarse_lattices is a list of transmission power
# fine_lattices is a list of power of fine lattices
# R is the rate constraint of the second hop (a list)
# A is a coefficient matrix defined on integer ring
# mod_scheme == 'opt_mod', or 'sim_mod'
def second_hop(trans_fine_lattices, trans_coarse_lattices, A, rate_sec_hop, mod_scheme):
    (M, L) = (A.nrows(), A.ncols())
    if M != L:
        raise Exception("L and M should be the same in destination's perspective.")
    # determine the actual finest lattice in m-th relay's linear combination
    relay_fine_lattices = [float('inf')]*M
    for i_M in range(0, M):
        for i_L in range(0, L):
            if (A[i_M, i_L] != 0) and (trans_fine_lattices[i_L]<relay_fine_lattices[i_M]):
                relay_fine_lattices[i_M] = trans_fine_lattices[i_L]
    
    # determine which lattices the m-th relay can do
    # modulo operation(constrained by forwarding rate)
    list_mod = []
    for i_M in range(0, M):
        list_mod_i_M = []
        for i_L in range(0, L):
            if rate_sec_hop[i_M] >= 0.5*log(trans_coarse_lattices[i_L]/relay_fine_lattices[i_M], 2):
                list_mod_i_M += [trans_coarse_lattices[i_L]]
        list_mod += [list_mod_i_M]
    
    if mod_scheme == 'opt_mod':
        # iterate all permutations A^L_L. If a permutation cannot be generated
        # with the lists above, then discard it.
        for mod_order in itertools.permutations(list(range(0, M)), L):
            # each element in mod_order: the relay that the l-th lattice should be assigned to
            mod_order_list = list(mod_order)
            is_valid = True
            for i_L in range(0, L):
                if trans_coarse_lattices[i_L] not in list_mod[mod_order_list[i_L]]:
                    # if not in, set as invalid
                    is_valid = False
                    break
            if is_valid == True:
                # if valid, test decodability
                is_decodable = True
                A_F = matrix(GF(p), A)
                A_i_L = A_F
                x = list(trans_coarse_lattices) # copy
                # del_col = []
                # del_row = []
                for i_L in range(0, L):
                    if A_i_L.rank() < L-i_L:
                        is_decodable = False
                        break
                    elif A_i_L.rank() == L-i_L:
                        x_min_idx, x_min_val = min(enumerate(x), key=lambda x:x[1])
                        x.pop(x_min_idx)
                        if A_i_L.rank() > 1:
                            # del_col += [x_min_idx]
                            # del_row += [mod_order[x_min_idx]]
                            # A_i_L = A_F.delete_columns(del_col)
                            # A_i_L = A_i_L.delete_rows(del_row)
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
                    return True # success
            else:
                # if invalid, then discard
                pass
        # if no successful scheme found, return failure
        return False # failure
    elif mod_scheme == 'sim_mod':
        # determine the coarsest lattice in each relay
        relay_coarsest_lattices = [float(0)]*M
        for i_M in range(0, M):
            for i_L in range(0, L):
                if (A[i_M, i_L] != 0) and (trans_coarse_lattices[i_L]>relay_coarsest_lattices[i_M]):
                    relay_coarsest_lattices[i_M] = trans_coarse_lattices[i_L]
        # if the m-th relay cannot mod relay_coarsest_lattices[i_M], then failure occurs.
        for i_M in range(0, M):
            if relay_coarsest_lattices[i_M] not in list_mod[i_M]:
                return False # failure
        # if all relay can forward, then the system can forward successfully
        return True # success 
    elif mod_scheme == 'naive_mod':
        coarsest_lattices = max(trans_coarse_lattices)
        # if the m-th relay cannot mod coarsest_lattices, then failure occurs.
        for i_M in range(0, M):
            if coarsest_lattices not in list_mod[i_M]:
                return False # failure
        return True
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
    fine_lattices = [0.5, 0.2]
    coarse_lattices = [1.5, 2.1]
    R = [1.5, 1.7]
    A = matrix(ZZ, 2, 2, [[1, 2], [1, 1]])
    p = 3 
    print second_hop(fine_lattices, coarse_lattices, A, R, 'opt_mod')
    
    fine_lattices = [0.5, 0.2, 0.4]
    coarse_lattices = [1.5, 2.1, 1]
    R = [1.5, 1.7, 3]
    A = matrix(ZZ, 3, 3, [[0, 0, 2], [2, 2, 0], [0, 2, 2]])
    print second_hop(fine_lattices, coarse_lattices, A, R, 'opt_mod')
    
    