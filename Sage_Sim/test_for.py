'''
Simulation of Compute and Forward with arbitrary coarse and fine
lattices.
Author: Yihua Tan
Email: adherentx.tyh@gmail.com
The Chinese University of Hong Kong
'''
import copy
import sys
from sage.all import *

import time
import itertools

print 'Hello, I am testing rank failure in Sage.'

N = 3
M = 3
L = 3
D = 2

for mod_order in itertools.permutations(list(range(0, M)), L):
    print mod_order



'''
for i_s in range(0, D**M):
            for i_a in range(0, M):
                print (i_s%(D**(i_a+1)))/(D**(i_a))
'''



'''
iter = 50000
rank_failure_list = [];
P = Primes()
p = P.first()
while true:
    if p > 100: break
    print p
    rank_failure = 0
    ring_c = IntegerModRing(p)
    for i in range(0, iter):
        A = matrix.random(ring_c, N, N)
        if A.rank() != N:
            rank_failure += 1
    rank_failure_list += [rank_failure]
    p = P.next(p)
# for i = range(0, iter)
rf_plot = list_plot(zip(primes_first_n(len(rank_failure_list)), [float(x)/iter for x in rank_failure_list]), \
                    axes_labels=['prime number $p$', 'failure probability'], \
                    title='the rank failure of a $2 \\times 2$ random matrix defined over GF(p)')
rf_plot.show()
rf_plot.save('/home/adherentx/Dropbox/Research/My_Report/Compute_and_Forward/Sage_Sim/rank_failure.eps')
raw_input()
'''