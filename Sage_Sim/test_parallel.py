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

print 'Hello, I am testing parallel computing in Sage.'

L = 2 # L transmitters
M = 2 # M relays
p = 7 # The prime number

@parallel(ncpus=3)
def par_func(a, b):
    set_random_seed()
    print matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
    return a*b

t1 = time.ctime()

a = [1, 6, 3]
b = [-1, 1, 1]
fac_results = sorted(list(par_func(zip(a, b))))
print fac_results
#print fac_results

t2 = time.ctime()

print 'Simulation started at ', t1
print 'Simulation ended at  ', t2

