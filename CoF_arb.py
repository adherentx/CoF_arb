'''
Simulation of Compute and Forward. In alpha test. 
Author: Yihua Tan
Email: adherentx.tyh@gmail.com
The Chinese University of Hong Kong
'''

import copy
import sys
from sage.all import *
from numpy import arange
import time
from CoF_LLL import *
from CoF_basic import *
from CoF_second_hop import *
from alternate_opt import *
from scipy import optimize
import math
import pickle
#import gc
#sys.path.insert(1, '/usr/local/lib/python2.7/dist-packages/')
#from guppy import hpy
#my_hpy = hpy()

print 'Hello, this is the simulation of CoF.'

# This is for test.

def CoF_compute_search_pow_flex(P_con, H_a, is_dual_hop, rate_sec_hop=[], mod_scheme='sym_mod', quan_scheme='sym_quan', beta=[]):
    (M, L) = (H_a.nrows(), H_a.ncols())
    if beta == []:
        beta = vector(RR, [1]*L)
    cof_pow = lambda x: -CoF_compute_fixed_pow_flex(x, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, beta)
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
            sum_rate_opt = CoF_compute_fixed_pow_flex(P_opt, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, beta)
        elif P_Search_Alg == 'anneal':
            res_cof = optimize.anneal(cof_pow, initial_guess, schedule='cauchy', T0=1, Tf=1e-6, \
                      full_output=True, maxiter=30, lower=[1, 1], upper=[P_con, P_con], dwell=30, disp=True)
            P_opt = list(res_cof[0])
            sum_rate_opt = -res_cof[1]
        elif P_Search_Alg == 'brute_fmin':
            res_brute = optimize.brute(cof_pow, Pranges, Ns=brute_fmin_number, full_output=True, finish=None)
            P_brute_opt = res_brute[0]
            sum_rate_brute = -res_brute[1] # negative! see minus sign in cof_pow
            res_fmin = optimize.fmin(cof_pow, P_brute_opt, xtol=0.1, ftol=0.001, maxiter=brute_fmin_maxiter, full_output=True)
            P_fmin_opt = res_fmin[0]
            sum_rate_opt = -res_fmin[1]
        else:
            raise Exception('error: algorithm not supported')
    except:
        print 'error in search algorithms'
        raise
    return sum_rate_opt
    

@parallel(ncpus=Cores)
def CoF_compute_eq_pow_con_first_hop(P_con, M, L):
    sum_rate = 0
    sum_rate_var = 0
    for i_H in range(0, batch_H):
        set_random_seed() # to avoid producing the same H_a in different threads
        if DEBUG_H == True:
            H_a = matrix(RR, M, L, [[0.430020183081381, 0.507485287402749], [0.571045084182643, -0.846256047736565]])
        else:
            H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
        
        sum_rate_i_H = 0
        sum_rate_i_H_var = 0
        
        # Fixed power
        sum_rate_A = CoF_compute_fixed_pow_flex((P_con, )*L, False, H_a, False)
        
        # Variable power
        if is_alternate == True:
            sum_rate_A_var = alternate_optimize(P_con, H_a, is_dual_hop=False)
        else:
            sum_rate_A_var = CoF_compute_search_pow_flex(P_con, H_a, is_dual_hop=False)

        sum_rate_i_H = sum_rate_A
        sum_rate_i_H_var = sum_rate_A_var
        
        sum_rate += sum_rate_i_H
        sum_rate_var += sum_rate_i_H_var
    # for i_H
    sum_rate /= batch_H
    sum_rate_var /= batch_H
    
#     print P_con, sum_rate, sum_rate_var
    return {'sum_rate': sum_rate, sum_rate_var: 'sum_rate_var'}

@parallel(ncpus=Cores)
def CoF_compute_eq_pow_con_dual_hops(P_con, M, L):
    sum_rate_fixed_pow_sym_mod = 0
    sum_rate_sym_mod = 0
    sum_rate_asym_mod = 0
    sum_rate_asym_mod_asym_quan = 0
    '''How to determine the power of the relays?'''
    P_relay = 0.25*P_con
    #P_relay = P_con
    #P_relay = sqrt(P_con) 
    for i_H in range(0, batch_H):
        if is_set_H == True:
            H_a = set_H_a
            H_b = set_H_b
        else:
            set_random_seed() # to avoid producing the same H_a in different threads
            H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
            set_random_seed()
            H_b = (matrix.random(RR, M, 1, distribution=RealDistribution('gaussian', 1))).column(0)
        if is_set_beta == True:
            beta = set_beta
        else:
            beta = vector(RR, [1]*L)
        rate_sec_hop = [0]*M # ? bit/s for each parallel channel in the second hop
        for i_h_b in range(0, M):
            rate_sec_hop[i_h_b] = 0.5*log(1+H_b[i_h_b]**2*P_relay, 2)
        
        
        try:
            (sum_rate_i_H_fixed_pow_sym_mod, sum_rate_i_H_sym_mod, sum_rate_i_H_asym_mod, sum_rate_i_H_asym_mod_asym_quan) = \
                CoF_compute_eq_pow_con_dual_hops_fixed_H(is_alternate, P_con, H_a, True, rate_sec_hop, beta)
        except:
            print 'error in searching for good power'
            raise
        
        sum_rate_fixed_pow_sym_mod += sum_rate_i_H_fixed_pow_sym_mod
        sum_rate_sym_mod += sum_rate_i_H_sym_mod
        sum_rate_asym_mod += sum_rate_i_H_asym_mod
        sum_rate_asym_mod_asym_quan += sum_rate_i_H_asym_mod_asym_quan
    # for i_H
    sum_rate_fixed_pow_sym_mod /= batch_H
    sum_rate_sym_mod /= batch_H
    sum_rate_asym_mod /= batch_H
    sum_rate_asym_mod_asym_quan /= batch_H
    
    return {'sum_rate_fixed_pow_sym_mod': sum_rate_fixed_pow_sym_mod,
            'sum_rate_sym_mod': sum_rate_sym_mod, 
            'sum_rate_asym_mod': sum_rate_asym_mod, 
            'sum_rate_asym_mod_asym_quan': sum_rate_asym_mod_asym_quan}


@parallel(ncpus=Cores)
def CoF_beta_search_dual_hops(P_con, M, L, beta):
    H_a = set_H_a
    H_b = set_H_b
    rate_sec_hop = [0]*M # ? bit/s for each parallel channel in the second hop
    P_relay = 0.25*P_con
    for i_h_b in range(0, M):
        rate_sec_hop[i_h_b] = 0.5*log(1+H_b[i_h_b]**2*P_relay, 2)
    
    (sum_rate_i_H_fixed_pow_sym_mod, sum_rate_i_H_sym_mod, sum_rate_i_H_asym_mod, sum_rate_i_H_asym_mod_asym_quan) = \
        CoF_compute_eq_pow_con_dual_hops_fixed_H(is_alternate, P_con, H_a, True, rate_sec_hop, beta)
    
    return {'sum_rate_i_H_fixed_pow_sym_mod': sum_rate_i_H_fixed_pow_sym_mod, 
            'sum_rate_i_H_sym_mod': sum_rate_i_H_sym_mod, 
            'sum_rate_i_H_asym_mod': sum_rate_i_H_asym_mod, 
            'sum_rate_i_H_asym_mod_asym_quan': sum_rate_i_H_asym_mod_asym_quan}
    

def CoF_compute_eq_pow_con_dual_hops_fixed_H(is_alternate, P_con, H_a, is_dual_hop, rate_sec_hop, beta):
    # Fixed power
    sum_rate_i_H_fixed_pow_sym_mod = CoF_compute_fixed_pow_flex_fine_lattice((P_con, )*L, H_a, rate_sec_hop, beta)
    
    # Variable power
    if is_alternate == True:
        sum_rate_i_H_sym_mod = alternate_optimize(P_con, H_a, True, rate_sec_hop, 'sym_mod', 'sym_quan', beta)
        sum_rate_i_H_asym_mod = alternate_optimize(P_con, H_a, True, rate_sec_hop, 'asym_mod', 'sym_quan', beta)
        sum_rate_i_H_asym_mod_asym_quan = 0 # TODO
    else:
        sum_rate_i_H_sym_mod = CoF_compute_search_pow_flex(P_con, H_a, True, rate_sec_hop, 'sym_mod', 'sym_quan', beta)
        sum_rate_i_H_asym_mod = CoF_compute_search_pow_flex(P_con, H_a, True, rate_sec_hop, 'asym_mod', 'sym_quan', beta)
        sum_rate_i_H_asym_mod_asym_quan = CoF_compute_search_pow_flex(P_con, H_a, True, rate_sec_hop, 'asym_mod', 'asym_quan', beta)
    return (sum_rate_i_H_fixed_pow_sym_mod, sum_rate_i_H_sym_mod, sum_rate_i_H_asym_mod, sum_rate_i_H_asym_mod_asym_quan)



if __name__ == "__main__": 
    '''Equal Power Constraint'''
    P_eq_dB_Min = float(P_MIN)
    P_eq_dB_Max = float(P_MAX)
    P_delta = 5
    P_eq_dB = arange(P_eq_dB_Min, P_eq_dB_Max, P_delta)
    P_eq = [10**(P_eq_dB_i/10) for P_eq_dB_i in P_eq_dB]
    Pl_con = P_eq
    
    num_batch = iter_H/batch_H # number of batch
    assert mod(num_batch, Cores)==0, 'For the sake of speed, the number of batches should be a multiple of Cores.'
    
    t1 = time.ctime()
    if is_alternate == True:
        dir_name = '/home/adherentx/Dropbox/Research/My_Report/CoF_arb/Sage_Sim/Simulation_Results/'+t1+'-'+'alternate'+'-M=L='+str(M)+'-iter_H='+str(iter_H)
    else:
        dir_name = '/home/adherentx/Dropbox/Research/My_Report/CoF_arb/Sage_Sim/Simulation_Results/'+t1+'-'+P_Search_Alg+'-M=L='+str(M)+'-iter_H='+str(iter_H)
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    original_dir = os.getcwd()
    os.chdir(dir_name) # change to the directory where simulation results should be placed
    
    '''First Hop'''
    if False:
        sum_rate = [0]*len(Pl_con)
        sum_rate_var = [0]*len(Pl_con)
        for i_P in range(0, len(Pl_con)):
            print 'First Hop Simulation: P_con=', Pl_con[i_P]
            result_bat = list(CoF_compute_eq_pow_con_first_hop([(Pl_con[i_P], M, L)]*num_batch))
            sum_rate_bat = [result_bat[i][1]['sum_rate'] for i in range(0, num_batch)]
            sum_rate_var_bat = [result_bat[i][1]['sum_rate_var'] for i in range(0, num_batch)]
            sum_rate[i_P] = sum(sum_rate_bat)/num_batch
            sum_rate_var[i_P] = sum(sum_rate_var_bat)/num_batch
            
        sum_rate = [RR(sum_rate[i]) for i in range(0, len(Pl_con))]
        sum_rate_var = [RR(sum_rate_var[i]) for i in range(0, len(Pl_con))]
        t2 = time.ctime()
        
        print 'Simulation of the first hop started at ', t1
        print 'Simulation of the first hop ended at ', t2
                
        plot_sum_rate = list_plot(zip(P_eq_dB, sum_rate), plotjoined=True, marker='o', \
                                  rgbcolor=Color('red'), linestyle="--", \
                                  legend_label= 'Sum rate with fixed power', \
                                  title = 'Comparison of Fixed and Variable Power Methods in the First Hop')
        plot_sum_rate_var = list_plot(zip(P_eq_dB, sum_rate_var), plotjoined=True, marker='x', \
                                      rgbcolor=Color('blue'), linestyle='-', \
                                      legend_label = 'Sum rate with variable power')
        plot_compare = plot_sum_rate+plot_sum_rate_var
        plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
        plot_compare.save('Comparison_Fixed_and_Variable_Power_in_the_First_Hop-' \
                          +P_Search_Alg+'-is_alternate='+str(is_alternate)+'-M=L='+str(M)+'.eps')
        show(plot_compare)
        pickle.dump((P_eq_dB, {'sum_rate': sum_rate, 
                               'sum_rate_var': sum_rate_var}), 
                    open('First_Hop.pkl', 'w'))
    
    '''Dual Hops'''
    if True:
        sum_rate_fixed_pow_sym_mod = [0]*len(Pl_con)
        sum_rate_sym_mod = [0]*len(Pl_con)
        sum_rate_asym_mod = [0]*len(Pl_con)
        sum_rate_asym_mod_asym_quan = [0]*len(Pl_con)
        t3 = time.ctime()
        # debug:
        #CoF_compute_eq_pow_con_dual_hops(100)
        
        for i_P in range(0, len(Pl_con)):
            print 'Dual Hop Simulation: P_con=', Pl_con[i_P]
            result_bat = list(CoF_compute_eq_pow_con_dual_hops([(Pl_con[i_P], M, L)]*num_batch))
            
            sum_rate_fixed_pow_sym_mod_bat = [result_bat[i][1]['sum_rate_fixed_pow_sym_mod'] for i in range(0, num_batch)]
            sum_rate_sym_mod_bat = [result_bat[i][1]['sum_rate_sym_mod'] for i in range(0, num_batch)]
            sum_rate_asym_mod_bat = [result_bat[i][1]['sum_rate_asym_mod'] for i in range(0, num_batch)]
            sum_rate_asym_mod_asym_quan_bat = [result_bat[i][1]['sum_rate_asym_mod_asym_quan'] for i in range(0, num_batch)]
            
            sum_rate_fixed_pow_sym_mod[i_P] = sum(sum_rate_fixed_pow_sym_mod_bat)/num_batch
            sum_rate_sym_mod[i_P] = sum(sum_rate_sym_mod_bat)/num_batch
            sum_rate_asym_mod[i_P] = sum(sum_rate_asym_mod_bat)/num_batch
            sum_rate_asym_mod_asym_quan[i_P] = sum(sum_rate_asym_mod_asym_quan_bat)/num_batch
            
        sum_rate_fixed_pow_sym_mod = [RR(sum_rate_fixed_pow_sym_mod[i]) for i in range(0, len(Pl_con))]
        sum_rate_sym_mod = [RR(sum_rate_sym_mod[i]) for i in range(0, len(Pl_con))]
        sum_rate_asym_mod = [RR(sum_rate_asym_mod[i]) for i in range(0, len(Pl_con))]
        sum_rate_asym_mod_asym_quan = [RR(sum_rate_asym_mod_asym_quan[i]) for i in range(0, len(Pl_con))]
        t4 = time.ctime()
        
        print 'Simulation of dual hops started at ', t3
        print 'Simulation of dual hops ended at ', t4
        
        plot_sum_rate_fixed_pow_sym_mod = list_plot(zip(P_eq_dB, sum_rate_fixed_pow_sym_mod), plotjoined=True, marker='D', \
                                  rgbcolor=Color('green'), linestyle="-.", \
                                  legend_label= 'Fixed power(Original CoF)', \
                                  title = 'Comparison in Two-Hop System')
        plot_sum_rate_sym_mod = list_plot(zip(P_eq_dB, sum_rate_sym_mod), plotjoined=True, marker='H', \
                                  rgbcolor=Color('black'), linestyle=":", \
                                  legend_label= 'Symmetric modulo approach with variable power')
        plot_sum_rate_asym_mod = list_plot(zip(P_eq_dB, sum_rate_asym_mod), plotjoined=True, marker='x', \
                                      rgbcolor=Color('blue'), linestyle='-', \
                                      legend_label = 'Asymmetric modulo approach with variable power')
        plot_sum_rate_asym_mod_asym_quan = list_plot(zip(P_eq_dB, sum_rate_asym_mod_asym_quan), plotjoined=True, marker='o', \
                                        rgbcolor=Color('red'), linestyle='--', \
                                        legend_label = 'Asymmetric modulo and asymmetric quantization approach with variable power')
        
        plot_compare = plot_sum_rate_fixed_pow_sym_mod+plot_sum_rate_sym_mod+plot_sum_rate_asym_mod+plot_sum_rate_asym_mod_asym_quan
        plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
        plot_compare.set_legend_options(loc='upper left')
        plot_compare.save('Comparison_in_Dual_Hops_System-' \
                          +P_Search_Alg+'-is_alternate='+str(is_alternate)+'-M=L='+str(M)+'.eps')
        show(plot_compare)
        pickle.dump((P_eq_dB, {'sum_rate_fixed_pow_sym_mod': sum_rate_fixed_pow_sym_mod, 
                               'sum_rate_sym_mod': sum_rate_sym_mod, 
                               'sum_rate_asym_mod': sum_rate_asym_mod, 
                               'sum_rate_asym_mod_asym_quan': sum_rate_asym_mod_asym_quan}), 
                    open('Dual_Hops.pkl', 'w'))
        
        print 'fixed power: '; print sum_rate_fixed_pow_sym_mod
        print 'symmetric mod with variable power: '; print sum_rate_sym_mod
        print 'asymmetric mod with variable power: '; print sum_rate_asym_mod
        print 'asymmetric mod and asymmetric quantization with variable power: '; print sum_rate_asym_mod_asym_quan
        


    '''Beta optimization'''
    beta_min = 0.2
    beta_max = 2
    beta_delta = 0.2
    beta_number = int(round((beta_max-beta_min)/beta_delta)+1)
    #assert beta_number%Cores==0, 'wrong beta_number setting!'
    P_con = 10**(20/10) # 40 dB
    result = []
    if False:
        if L==M==2:
            result = list(CoF_beta_search_dual_hops([(P_con, M, L, vector(RR, [beta_min+j0*beta_delta, beta_min+j1*beta_delta])) for j0 in range(0,beta_number) for j1 in range(0,beta_number)]))
        else:
            pass
        # (sum_rate_i_H_fixed_pow_sym_mod, sum_rate_i_H_sym_mod, sum_rate_i_H_asym_mod, sum_rate_i_H_asym_mod_asym_quan)
        for re in result:
            print re
        pickle.dump(result, open('beta_search.pkl', 'w'))
        best_result_0 = max(result, key=lambda x: x[1]['sum_rate_i_H_asym_mod_asym_quan'])
        print 'The best result for asym_mod_asym_quan is: \n', best_result_0
        best_result_1 = max(result, key=lambda x: x[1]['sum_rate_i_H_asym_mod'])
        print 'The best result for asym_mod is: \n', best_result_1
        best_result_2 = max(result, key=lambda x: x[1]['sum_rate_i_H_sym_mod'])
        print 'The best result for sym_mod is: \n', best_result_2
        best_result_3 = max(result, key=lambda x: x[1]['sum_rate_i_H_fixed_pow_sym_mod'])
        print 'The best result for fixed_pow_sym_mod is: \n', best_result_3
        pickle.dump({'sum_rate_i_H_asym_mod_asym_quan': best_result_0, 
                     'sum_rate_i_H_asym_mod': best_result_1, 
                     'sum_rate_i_H_sym_mod': best_result_2, 
                     'sum_rate_i_H_fixed_pow_sym_mod': best_result_3},
                    open('best_beta.pkl', 'w'))
    # P_con, M, L, beta
    
    os.chdir(original_dir) # recover directory
    raw_input() # stop Sage from shutting down
