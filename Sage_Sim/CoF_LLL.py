'''
This is something related with lattice reduction algorithms.
'''
from sage.all import *
import copy
from numpy import arange
from sage.parallel.all import *
import time
from CoF_basic import *
from CoF_second_hop import *


def CoF_compute_fixed_pow_flex_fine_lattice(P_t, H_a, rate_sec_hop):
    M = 2
    L = 2
    P1, P2 = P_t
    if math.isnan(P1) or math.isnan(P2):
        print 'P1 or P2 should not be NaN!'
        return 0
    if P1 <= 0 or P2 <= 0:
        print 'P1 and P2 should be positive'
        return 0
    P_vec = vector(RR, [P1, P2])
    P_mat = matrix.diagonal([sqrt(x) for x in P_vec])
    # Use LLL to find a good A matrix
    # determine the fine lattice of m-th relay at the same time
    try:
        (A_best_LLL, rate_list_A_LLL, relay_fine_lattices) = Find_A_and_Rate(P_mat, P_vec, H_a, is_return_rate_list=True)
    except:
        print 'error in seeking A and rate'
        raise
    rate_list = list(rate_list_A_LLL)
    for i_l in range(0, L):
        for i_m in range(0, M):
            if (A_best_LLL[i_m, i_l]!=0) and (rate_sec_hop[i_m]<rate_list[i_l]):
                rate_list[i_l] = rate_sec_hop[i_m]
    sum_rate = sum(rate_list)
    return sum_rate

'''This is for L=M=2!'''
def CoF_compute_fixed_pow(P_t, is_return_A, *params):
    P1, P2 = P_t
    #print 'P1 = ', P1, '   P2 = ', P2
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
        # return (A_best_LLL, sum_rate_A_LLL, alpha_opt_LLL)
        if is_return_A == True:
            return (sum_rate_A_LLL, A_best_LLL)
        else:
            return sum_rate_A_LLL
    else:
        # return (zero_matrix(ZZ, M, L), 0, zero_vector(RR, M))
        if is_return_A == True:
            return (0, zero_matrix(ZZ, M, L))
        else:
            return 0
        

def Find_A_m_list(P, H):
    # P is a LxL matrix, H is a MxL matrix
    (M, L) = (H.nrows(), H.ncols())
    A_m_list = []
    for i_a in range(0, M):
        h = H.row(i_a)
        G = P*(identity_matrix(RR, L)-(P.T*h.column()*h.row()*P/(1+((h.row()*P).norm(p=2))**2)))*P.T
        G_rdf = G.change_ring(RDF)
        if G_rdf.is_positive_definite() == False:
            raise Exception('G should be positive definite!')
        F = G_rdf.cholesky()
        if F.rank() != L:
            raise Exception('F should be full rank. Check it!')
        amp = 10**5
        F_a = amp*F
        F_a_rdf = F_a.change_ring(RDF)
        F_a_z = matrix(ZZ, L, L, [round(x) for x in F_a_rdf.list()])
        F_a_reduced = F_a_z.LLL(use_givens=True)
        F_reduced = F_a_reduced/amp
        F_reduced_list = F_reduced.rows()
        F_reduced_list = sorted(F_reduced_list, key=lambda x:x.norm(p=2), reverse=False)
        F_reduced = matrix(F_reduced_list)
        FF = F_reduced*F.inverse()
        FF_z = matrix(ZZ, L, L, [round(x) for x in FF.list()])
        A_m_list += [FF_z]
    return A_m_list

def Find_A_and_Rate(P_mat, P_vec, H, is_return_rate_list=False):
    # Use LLL to find a good A matrix
    (M, L) = (H.nrows(), H.ncols())
    A_list = Find_A_m_list(P_mat, H)
    A = zero_matrix(ZZ, M, L)
    for i_a in range(0, M):
        A.set_row(i_a, A_list[i_a].row(0))
    A_F = matrix(GF(p), A)
    rank_first_row = A_F.rank()
#     alpha_opt = zero_vector(RR, M)
    relay_fine_lattices = [0]*M
    sum_rate = 0
    rate_list = [0]*L
    A_best = zero_matrix(ZZ, M, L)
    if rank_first_row == min(L, M):
        # full rank
#         for i_alpha in range(0, M):
#             alpha_opt[i_alpha] = alpha_find(H.row(i_alpha), P_mat, A.row(i_alpha))
#         rate = rate_computation(L, M, P_vec, alpha_opt, H, A)
        try:
            (rate, relay_fine_lattices) = rate_computation_MMSE_alpha(L, M, P_vec, H, A)
        except:
            print 'error in rate_computation_MMSE_alpha: pos 0'
            raise
        sum_rate = sum(rate)
        rate_list = list(rate)
        A_best = A
    else:
        # rank deficiency: search for full-rank matrix A
        D = A.nrows() # search in [0, D-1] in each matrix in the list
        for i_s in range(0, D**M):
            A = zero_matrix(RR, M, L)
            for i_a in range(0, M):
                idx_row_i_a = (i_s%(D**(i_a+1)))/(D**i_a)
                A.set_row(i_a, A_list[i_a].row(idx_row_i_a))
            A_F = matrix(GF(p), A)
            if A_F.rank() == min(L, M):
                try:
                    (rate, relay_fine_lattices_i_row_search) = rate_computation_MMSE_alpha(L, M, P_vec, H, A)
                except:
                    print 'error in rate_computation_MMSE_alpha: pos 1'
                if sum_rate < sum(rate):
                    sum_rate = sum(rate)
                    rate_list = list(rate)
#                     alpha_opt = alpha_opt_i_row_search
                    relay_fine_lattices = relay_fine_lattices_i_row_search
                    A_best = A
            else:
                pass
#     A_best_F = matrix(GF(p), A_best)
#     if A_best_F.rank() < min(L, M):
#         pass # for test
#     if A_best == matrix(GF(p), 2, 2, [[1, 2], [2, 1]]):
#         pass
    if is_return_rate_list == True:
        return (A_best, rate_list, relay_fine_lattices)
    else:
        return (A_best, sum_rate, relay_fine_lattices)

    
@parallel(ncpus=Cores)
def CoF_rank_deficiency(P_con):
    iter_H = 2000
    rank_deficiency = 0
    for i_H in range(0, iter_H):
        set_random_seed() # to avoid producing the same H in different threads
        H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
        # H_a = matrix(RR, M, L, [[-0.869243498832414, 0.762538785616652], [-0.914996809982352, 0.801032403084523]])
        # Fixed power
        P_mat = sqrt(P_con)*identity_matrix(L)
        P_vec = P_con*vector([1]*L)
        
        # Use LLL to find a good A matrix
        (A, dummy, dummy) = Find_A_and_Rate(P_mat, P_vec, H_a, is_return_rate_list=False)
        
        A_F = matrix(GF(p), A)
        alpha_opt = zero_vector(RR, M)
        if A_F.rank() == min(L, M):
            pass
        else:
            # rank deficiency: outage
            sum_rate_A = 0
            rank_deficiency += 1
    
    return float(rank_deficiency)/iter_H


if __name__ == "__main__":
    P_eq_dB_Min = float(20)
    P_eq_dB_Max = float(40)
    P_delta = 5
    P_eq_dB = arange(P_eq_dB_Min, P_eq_dB_Max, P_delta)
    P_eq = [10**(P_eq_dB_i/10) for P_eq_dB_i in P_eq_dB]
    Pl_con = P_eq
    
    t1 = time.ctime()
    result = list(CoF_rank_deficiency(Pl_con))
    t2 = time.ctime()
    
    print 'Simulation started at ', t1
    print 'Simulation ended at ', t2
    rank_deficiency = [result[i][1] for i in range(0, len(Pl_con))]
    plot_rank_deficiency = list_plot(zip(P_eq_dB, rank_deficiency), plotjoined=True, marker='o', \
                              rgbcolor=Color('red'), linestyle="--", \
                              legend_label= 'probability', \
                              title = 'Rank deficiency probability')
    print 'rank_deficiency = ', rank_deficiency
    show(plot_rank_deficiency)
    
    raw_input() # stop Sage from shutting down
    
    
