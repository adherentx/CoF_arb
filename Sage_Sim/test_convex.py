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

print 'Hello, I am testing concavity.'

P = matrix(SR, 2, 2, [var('P1'), 0, 0, var('P2')])
A = matrix(SR, 2, 2, var('a11', 'a12', 'a21', 'a22'))
H = matrix(SR, 2, 2, var('h11', 'h12', 'h21', 'h22'))
#item1 = P1/(a11**2*P1+a12**2*P2-(h11*P1*a11+h12*P2*a12)**2/(1+h11**2*P1+h12**2*P2))
item1 = P1/(a21**2*P1+a22**2*P2-(h21*P1*a21+h22*P2*a22)**2/(1+h21**2*P1+h22**2*P2))
item2 = P2/(a21**2*P1+a22**2*P2-(h21*P1*a21+h22*P2*a22)**2/(1+h21**2*P1+h22**2*P2))

phi = log(item1*item2, 2)

# A1 = matrix.random(ZZ, 2, 2)
# H1 = matrix.random(RR, 2, 2)
A1 = matrix(ZZ, 2, 2, [[1, 1], [-1, 1]])
H1 = matrix(RR, 2, 2, [[0.430020183081381, 0.507485287402749], [0.571045084182643, -0.846256047736565]])
 
phi_p = phi(a11=A1[0,0], a12=A1[0,1], a21=A1[1,0], a22=A1[1,1], \
            h11=H1[0,0], h12=H1[0,1], h21=H1[1,0], h22=H1[1,1])
p = plot3d(phi_p, [P1, 0.01, 100], [P2, 0.01, 100], axes=True)
p.show()
p1 = plot(phi_p(P2=100), xmin=0.1, xmax=100)
p1.show()
raw_input()


