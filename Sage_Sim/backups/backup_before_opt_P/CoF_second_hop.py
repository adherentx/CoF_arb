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
            is_valid = True
            for i_L in range(0, L):
                if trans_coarse_lattices[i_L] not in list_mod[mod_order[i_L]]:
                    # if not in, set as invalid
                    is_valid = False
                    break
            if is_valid == True:
                # if valid, test decodability
                is_decodable = True
                A_i_L = matrix(GF(p), A)
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
                            A_i_L = A_i_L.delete_rows([mod_order[x_min_idx]])
                        else:
                            # don't need to reducing since no further operations need it
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
                
    else:
        raise Exception("mod_scheme should take value as 'opt_mod' or 'sim_opt'!")
    

            

if __name__ == "__main__":
    print '-----------------------------------\n'+ \
        'testing CoF_second_hop\n'
    fine_lattices = [0.5, 0.2]
    coarse_lattices = [1.5, 2.1]
    R = [1.5, 1.7]
    A = matrix(ZZ, 2, 2, [[1, 2], [1, 1]])
    print second_hop(fine_lattices, coarse_lattices, A, R)
    
    