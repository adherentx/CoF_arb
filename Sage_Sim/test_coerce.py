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
from scipy import optimize
print 'Hello, I am testing coercing.'

class my_class:
    def __init__(self, x, s):
        self.x = x
        self.s = s
    def __float__(self):
        return self.x

def f(x):
    return my_class(x**2-2*x+1, 'apple')

res_f_brute = optimize.brute(f, (slice(-4, 4, 1), ), full_output=True)

print res_f_brute