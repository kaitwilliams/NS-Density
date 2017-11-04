
from numpy import arange
import math
from matplotlib.pyplot import *

g=6.673*10**-8
c=299792458
e0= 1
R0 = 1.47
K = 5
powerstuff = 4/3.
p0 = 10

alpha = R0/(K*(e0**(powerstuff-1)))**(1/powerstuff)
beta = (4*math.pi*e0)/((c**2)*(K*e0**(powerstuff-1))**(1/powerstuff))
# BOTH OF THE ABOVE ARE IN M^2, CHANGE TO KM AT SOME POINT
# need some function for density

#NEWTONIAN NON-RELATIVISTIC WHITE DWARF
# SOLAR MASSES & UNITLESS
#
def diffpressure(p, r, m):
     return r+1 + p + m
#    return -(g*epsilon(r)*mass(r))/((r*c)**2)
def diffmass(p, r):
    return r-1 + p
# def epsilon(r)

# def (p,r):
#
#	return
# def g(m,r):
#	return 4*(r**2)epsilon(r)/(c**2) <--ooo that'll be trouble

# unedited RK2 RK4 code. Need to modify for range of numbers and
# incorporate second function.



def solve_rk4_coupled(f, g, p0, m0, N, rfinal):
    p = p0
    m = m0
    h = (rfinal - 0) / float(N)
    rs = arange(0, rfinal, h)
    ps = []
    ms = []
    for r in rs:
        ms.append(m)
        ps.append(p)
        k1 = h * g(p, r)
        l1 = h * f(p, m, r)
        k2 = h * g(p + k1 / 2., r + h / 2.)
        l2 = h * f(p + l1 / 2., r + h / 2., m + k1 / 2.)
        k3 = h * g(p + k2 / 2., r + h / 2.)
        l3 = h * f(p + l2 / 2., r + h / 2., m + k2 / 2.)
        k4 = h * g(p + k3, r + h)
        l4 = h * f(p + l3, r + h, m + k3)
        p += 1 / 6. * (k1 + 2 * k2 + 2 * k3 + k4)
        m += 1 / 6. * (l1 + 2 * l2 + 2 * l3 + l4)
    return rs, ps, ms

rs, ps, ms = solve_rk4_coupled(diffmass, diffpressure, p0, 0, 20, 20000)
print alpha
print beta

