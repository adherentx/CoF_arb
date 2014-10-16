from sage.all import *
sys.path.append('/home/adherentx/workspace/CoF_arb/Sage_Sim/')
from CoF_basic import *
import pickle
from numpy import arange

def result_combine_first_hop(result1, result2):
    sum_rate = result1.sum_rate + result2.sum_rate
    sum_rate_var = result1.sum_rate_var + result2.sum_rate_var
    return (sum_rate, sum_rate_var)


def result_combine_two_hop(result1, result2):
    sum_rate_fixed_pow_sim_mod = result1.sum_rate_fixed_pow_sim_mod + result2.sum_rate_fixed_pow_sim_mod
    sum_rate_sim_mod = result1.sum_rate_sim_mod + result2.sum_rate_sim_mod
    sum_rate_opt_mod = result1.sum_rate_opt_mod + result2.sum_rate_opt_mod
    return (sum_rate_fixed_pow_sim_mod, sum_rate_sim_mod, sum_rate_opt_mod)

if __name__ == '__main__':
    
    
    
    P_eq_dB_Min = float(10)
    P_eq_dB_Max = float(90)
    P_delta = 5
    P_eq_dB = arange(P_eq_dB_Min, P_eq_dB_Max, P_delta)
    P_eq = [10**(P_eq_dB_i/10) for P_eq_dB_i in P_eq_dB]
    Pl_con = P_eq
    P_Search_Alg = 'brute'
    is_alternate = 'False'
    
    # First Hop
    Result_First_Hop1 = pickle.load(open('First_Hop1.pkl', 'r'))
    Result_First_Hop2 = pickle.load(open('First_Hop2.pkl', 'r'))
    (sum_rate, sum_rate_var) = result_combine_first_hop(Result_First_Hop1, Result_First_Hop2)
    
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
    pickle.dump((P_eq_dB, CoF_Sim_Result(sum_rate, sum_rate_var)), open('First_Hop.pkl', 'w'))
    
    
    # Two hops
    Result_Two_Hop1 = pickle.load(open('Dual_Hops1.pkl', 'r'))
    Result_Two_Hop2 = pickle.load(open('Dual_Hops2.pkl', 'r'))
    (sum_rate_fixed_pow_sim_mod, sum_rate_sim_mod, sum_rate_opt_mod) = result_combine_two_hop(Result_Two_Hop1, Result_Two_Hop2)
    
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
    pickle.dump((P_eq_dB, CoF_Dual_Hops_Sim_Result(sum_rate_fixed_pow_sim_mod, sum_rate_sim_mod, sum_rate_opt_mod)), open('Dual_Hops.pkl', 'w'))
    
    
    