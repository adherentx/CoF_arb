'''
This is something related with lattice reduction algorithms.
'''
from sage.all import *
import copy
from numpy import arange
from sage.parallel.all import *
import time

L = 2 # L transmitters
M = 2 # M relays
p = 13 # The prime number

Cores = 1 # The number of CPU cores used in parallel computing
DEBUG_H = False # When this value is True, the channel matrix H_a is set as certain matrices


def Find_A(P, H):
    # P is a LxL matrix, H is a MxL matrix
    (M, L) = (H.nrows(), H.ncols())
    A = zero_matrix(ZZ, M, L)
    
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
        if F_a_reduced.rank() != 1:
            pass
            # raise Exception('The rank of F_a_reduced should be 1!')
        F_reduced = F_a_reduced/amp
        F_reduced_list = F_reduced.rows()
        F_reduced_list = sorted(F_reduced_list, key=lambda x:x.norm(), reverse=False)
        F_reduced = matrix(F_reduced_list)
        FF = F_reduced*F.inverse()
        FF_z = matrix(ZZ, L, L, [round(x) for x in FF.list()])
        idx = 0
        row_norm_idx = 0
        phi = float('inf')
        # In the LLL result, find the vector that can minimize phi
        for i_FF in range(0, L):
            phi_i_FF = (FF_z.row(i_FF)*G_rdf*FF_z.row(i_FF).column()).norm()
            if  phi_i_FF < phi:
                phi = phi_i_FF
                idx = i_FF
        if idx != 0:
            raise Exception('The best row should be the first one.')
        A.set_row(i_a, FF_z.row(idx))
    return A


@parallel(ncpus=Cores)
def CoF_rank_deficiency(P_con):
    iter_H = 2000
    rank_deficiency = 0
    for i_H in range(0, iter_H):
        set_random_seed() # to avoid producing the same H_a in different threads
        H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))

        # Fixed power
        P_mat = sqrt(P_con)*identity_matrix(L)
        P_vec = P_con*vector([1]*L)
        
        # Use LLL to find a good A matrix
        A = Find_A(P_mat, H_a)
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
    show(plot_rank_deficiency)
    
    raw_input() # stop Sage from shutting down
    
    
