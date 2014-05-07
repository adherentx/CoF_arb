'''
Simulation of Compute and Forward with arbitrary coarse and fine
lattices.
Author: Yihua Tan
Email: adherentx.tyh@gmail.com
The Chinese University of Hong Kong
'''


''' MEMO:
  Fix all the 'FIX ME' bugs!!!
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

def CoF_compute_search_pow(P_con, H_a, is_dual_hop, rate_sec_hop=[], mod_scheme='sim_mod'):
    cof_pow = lambda x: -CoF_compute_fixed_pow(x, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme)
    Pranges = ((0.1, P_con), (0.1, P_con))
    initial_guess = [0.5*P_con, 0.5*P_con]
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
    

@parallel(ncpus=Cores)
def CoF_compute_eq_pow_con_first_hop(P_con):
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
        sum_rate_A = CoF_compute_fixed_pow((P_con, P_con), False, H_a, False)
        
        # Variable power
        if is_alternate == True:
            sum_rate_A_var = alternate_optimize(P_con, H_a, is_dual_hop=False)
        else:
            sum_rate_A_var = CoF_compute_search_pow(P_con, H_a, is_dual_hop=False)

        sum_rate_i_H = sum_rate_A
        sum_rate_i_H_var = sum_rate_A_var
        
        sum_rate += sum_rate_i_H
        sum_rate_var += sum_rate_i_H_var
    # for i_H
    sum_rate /= batch_H
    sum_rate_var /= batch_H
    
#     print P_con, sum_rate, sum_rate_var
    return CoF_Sim_Result(sum_rate, sum_rate_var)

@parallel(ncpus=Cores)
def CoF_compute_eq_pow_con_dual_hops(P_con):
    sum_rate_fixed_pow_sim_mod = 0
    sum_rate_sim_mod = 0
    sum_rate_opt_mod = 0
    '''How to determine the power of the relays?'''
    P_relay = 0.25*P_con
    #P_relay = P_con
    #P_relay = sqrt(P_con) 
    for i_H in range(0, batch_H):
        set_random_seed() # to avoid producing the same H_a in different threads
        if DEBUG_H == True:
            H_a = matrix(RR, M, L, [[-0.333414283246484, 0.675839593022605], [0.000374794674703693, 0.766514412738790]])
        else:
            H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
            #H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 0.25))
        
        '''Fix Me: set the rate constraint of the second hop'''
        set_random_seed()
        H_b = (matrix.random(RR, M, 1, distribution=RealDistribution('gaussian', 1))).column(0)
        rate_sec_hop = [0]*M # [4, 6] # ? bit/s for each parallel channel in the second hop
        for i_h_b in range(0, M):
            rate_sec_hop[i_h_b] = 0.5*log(1+H_b[i_h_b]**2*P_relay, 2)
        
        # Fixed power
        sum_rate_i_H_fixed_pow_sim_mod = CoF_compute_fixed_pow_flex_fine_lattice((P_con, P_con), H_a, rate_sec_hop)
        
        # In dual-hop system, we use variable power method to avoid outage at the second hop.
        try:
            sum_rate_i_H_sim_mod = CoF_compute_search_pow(P_con, H_a, is_dual_hop=True, rate_sec_hop=rate_sec_hop, mod_scheme='sim_mod')
            sum_rate_i_H_opt_mod = CoF_compute_search_pow(P_con, H_a, is_dual_hop=True, rate_sec_hop=rate_sec_hop, mod_scheme='opt_mod')
        except:
            print 'error in searching for good power'
            raise
        
        sum_rate_fixed_pow_sim_mod += sum_rate_i_H_fixed_pow_sim_mod
        sum_rate_sim_mod += sum_rate_i_H_sim_mod
        sum_rate_opt_mod += sum_rate_i_H_opt_mod
    # for i_H
    sum_rate_fixed_pow_sim_mod /= batch_H
    sum_rate_sim_mod /= batch_H
    sum_rate_opt_mod /= batch_H
    
    return CoF_Dual_Hops_Sim_Result(sum_rate_fixed_pow_sim_mod, sum_rate_sim_mod, sum_rate_opt_mod)


if __name__ == "__main__": 
    '''Equal Power Constraint'''
    P_eq_dB_Min = float(10)
    P_eq_dB_Max = float(90)
    P_delta = 10
    P_eq_dB = arange(P_eq_dB_Min, P_eq_dB_Max, P_delta)
    P_eq = [10**(P_eq_dB_i/10) for P_eq_dB_i in P_eq_dB]
    Pl_con = P_eq
    
    num_batch = iter_H/batch_H # number of batch
    assert mod(num_batch, Cores)==0, 'For the sake of speed, the number of batches should be a multiple of Cores.'
    
    t1 = time.ctime()
    dir_name = '/home/adherentx/Dropbox/Research/My_Report/CoF_arb/Sage_Sim/Simulation_Results/'+t1+'-'+P_Search_Alg+'-iter_H='+str(iter_H)
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    original_dir = os.getcwd()
    os.chdir(dir_name) # change to the directory where simulation results should be placed
    
    '''First Hop'''
    if True:
        sum_rate = [0]*len(Pl_con)
        sum_rate_var = [0]*len(Pl_con)
        for i_P in range(0, len(Pl_con)):
            print 'First Hop Simulation: P_con=', Pl_con[i_P]
            result_bat = list(CoF_compute_eq_pow_con_first_hop([Pl_con[i_P]]*num_batch))
            sum_rate_bat = [result_bat[i][1].sum_rate for i in range(0, num_batch)]
            sum_rate_var_bat = [result_bat[i][1].sum_rate_var for i in range(0, num_batch)]
            sum_rate[i_P] = sum(sum_rate_bat)/num_batch
            sum_rate_var[i_P] = sum(sum_rate_var_bat)/num_batch
            
        sum_rate = [RR(sum_rate[i]) for i in range(0, len(Pl_con))]
        sum_rate_var = [RR(sum_rate_var[i]) for i in range(0, len(Pl_con))]
        t2 = time.ctime()
        
        print 'Simulation of the first hop started at ', t1
        print 'Simulation of the first hop ended at ', t2
                
        plot_sum_rate = list_plot(zip(P_eq_dB, sum_rate), plotjoined=True, marker='o', \
                                  rgbcolor=Color('red'), linestyle="--", \
                                  legend_label= 'sum rate of fixed power method', \
                                  title = 'Comparison of fixed and variable methods int the first hop')
        plot_sum_rate_var = list_plot(zip(P_eq_dB, sum_rate_var), plotjoined=True, marker='x', \
                                      rgbcolor=Color('blue'), linestyle='-', \
                                      legend_label = 'sum rate of variable power method')
        plot_compare = plot_sum_rate+plot_sum_rate_var
        plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
        plot_compare.save('Comparison_Fixed_and_Variable_Power_in_the_First_Hop-' \
                          +P_Search_Alg+'-is_alternate='+str(is_alternate)+'.eps')
        show(plot_compare)
        pickle.dump((P_eq_dB, CoF_Sim_Result(sum_rate, sum_rate_var)), open('First_Hop.pkl', 'w'))
    
    '''Dual Hops'''
    if False:
        sum_rate_fixed_pow_sim_mod = [0]*len(Pl_con)
        sum_rate_sim_mod = [0]*len(Pl_con)
        sum_rate_opt_mod = [0]*len(Pl_con)
        t3 = time.ctime()
        # debug:
        #CoF_compute_eq_pow_con_dual_hops(100)
        
        for i_P in range(0, len(Pl_con)):
            print 'Dual Hop Simulation: P_con=', Pl_con[i_P]
            result_bat = list(CoF_compute_eq_pow_con_dual_hops([Pl_con[i_P]]*num_batch))
            sum_rate_fixed_pow_sim_mod_bat = [result_bat[i][1].sum_rate_fixed_pow_sim_mod for i in range(0, num_batch)]
            sum_rate_sim_mod_bat = [result_bat[i][1].sum_rate_sim_mod for i in range(0, num_batch)]
            sum_rate_opt_mod_bat = [result_bat[i][1].sum_rate_opt_mod for i in range(0, num_batch)]
            sum_rate_fixed_pow_sim_mod[i_P] = sum(sum_rate_fixed_pow_sim_mod_bat)/num_batch
            sum_rate_sim_mod[i_P] = sum(sum_rate_sim_mod_bat)/num_batch
            sum_rate_opt_mod[i_P] = sum(sum_rate_opt_mod_bat)/num_batch
        
        sum_rate_fixed_pow_sim_mod = [RR(sum_rate_fixed_pow_sim_mod[i]) for i in range(0, len(Pl_con))]
        sum_rate_sim_mod = [RR(sum_rate_sim_mod[i]) for i in range(0, len(Pl_con))]
        sum_rate_opt_mod = [RR(sum_rate_opt_mod[i]) for i in range(0, len(Pl_con))]
        t4 = time.ctime()
        
        print 'Simulation of dual hops started at ', t3
        print 'Simulation of dual hops ended at ', t4
        
        plot_sum_rate_fixed_pow_sim_mod = list_plot(zip(P_eq_dB, sum_rate_fixed_pow_sim_mod), plotjoined=True, marker='D', \
                                  rgbcolor=Color('green'), linestyle="-.", \
                                  legend_label= 'simple modulo method with fixed power', \
                                  title = 'Comparison in dual-hops system')
        plot_sum_rate_sim_mod = list_plot(zip(P_eq_dB, sum_rate_sim_mod), plotjoined=True, marker='o', \
                                  rgbcolor=Color('red'), linestyle="--", \
                                  legend_label= 'simple modulo method with variable power')
        plot_sum_rate_opt_mod = list_plot(zip(P_eq_dB, sum_rate_opt_mod), plotjoined=True, marker='x', \
                                      rgbcolor=Color('blue'), linestyle='-', \
                                      legend_label = 'optimal mudulo method with variable power')
        plot_compare = plot_sum_rate_fixed_pow_sim_mod+plot_sum_rate_sim_mod+plot_sum_rate_opt_mod
        plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
        plot_compare.set_legend_options(loc='upper left')
        plot_compare.save('Comparison_in_Dual_Hops_System-' \
                          +P_Search_Alg+'-is_alternate='+str(is_alternate)+'.eps')
        show(plot_compare)
        pickle.dump((P_eq_dB, CoF_Dual_Hops_Sim_Result(sum_rate_fixed_pow_sim_mod, sum_rate_sim_mod, sum_rate_opt_mod)), open('Dual_Hops.pkl', 'w'))
    
    os.chdir(original_dir) # recover directory
    raw_input() # stop Sage from shutting down
