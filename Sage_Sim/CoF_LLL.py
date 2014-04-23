'''
This is something related with lattice reduction algorithms.
'''
from sage.all import *
import copy
from numpy import arange
from sage.parallel.all import *
import time
from CoF_basic import *


def alpha_find(h, P_mat, a):
    alpha_opt = (h.row()*P_mat*P_mat.T*a.column())[0,0]/(1+(h.row()*P_mat).norm()**2)
    return alpha_opt
    #alpha_opt_v = h.row()*P_mat*P_mat.T*a.column()*~(1+h.row()*P_mat*P_mat.T*h.column())
    #return alpha_opt_v[0, 0]

def Find_A_m_list(P, H):
    # P is a LxL matrix, H is a MxL matrix
    (M, L) = (H.nrows(), H.ncols())
    A_m_list = []
    for i_a in range(0, M):
        h = H.row(i_a)
        G = P*(identity_matrix(RR, L)-(P.T*h.column()*h.row()*P/(1+((h.row()*P).norm())**2)))*P.T
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
        F_reduced_list = sorted(F_reduced_list, key=lambda x:x.norm(), reverse=False)
        F_reduced = matrix(F_reduced_list)
        FF = F_reduced*F.inverse()
        FF_z = matrix(ZZ, L, L, [round(x) for x in FF.list()])
        A_m_list += [FF_z]
    return A_m_list

def Find_A_and_Rate(P_mat, P_vec, H):
    # Use LLL to find a good A matrix
    (M, L) = (H.nrows(), H.ncols())
    A_list = Find_A_m_list(P_mat, H)
    A = zero_matrix(ZZ, M, L)
    for i_a in range(0, M):
        A.set_row(i_a, A_list[i_a].row(0))
    A_F = matrix(GF(p), A)
    rank_first_row = A_F.rank()
    alpha_opt = zero_vector(RR, M)
    sum_rate = 0
    A_best = zero_matrix(ZZ, M, L)
    if rank_first_row == min(L, M):
        # full rank
        for i_alpha in range(0, M):
            alpha_opt[i_alpha] = alpha_find(H.row(i_alpha), P_mat, A.row(i_alpha))
        rate = rate_computation(L, M, P_vec, alpha_opt, H, A)
        sum_rate = sum(rate)
        A_best = A
    else:
        # rank deficiency: search for full-rank matrix A
        D = int(L-rank_first_row+1) # search in [0, D-1] in each matrix in the list
        for i_s in range(0, D**M):
            A = zero_matrix(RR, M, L)
            for i_a in range(0, M):
                idx_row_i_a = (i_s%(D**(i_a+1)))/(D**i_a)
                A.set_row(i_a, A_list[i_a].row(idx_row_i_a))
            A_F = matrix(GF(p), A)
            if A_F.rank() == min(L, M):
                alpha_opt_i_row_search = zero_vector(RR, M)
                for i_alpha in range(0, M):
                    alpha_opt_i_row_search[i_alpha] = alpha_find(H.row(i_alpha), P_mat, A.row(i_alpha))
                rate = rate_computation(L, M, P_vec, alpha_opt_i_row_search, H, A)
                if sum_rate < sum(rate):
                    sum_rate = sum(rate)
                    alpha_opt = alpha_opt_i_row_search
                    A_best = A
            else:
                pass
#     A_best_F = matrix(GF(p), A_best)
#     if A_best_F.rank() < min(L, M):
#         pass # for test
#     if A_best == matrix(GF(p), 2, 2, [[1, 2], [2, 1]]):
#         pass
    return (A_best, sum_rate, alpha_opt)

    
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
        (A, dummy, dummy) = Find_A_and_Rate(P_mat, P_vec, H_a)
        
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
    
    
