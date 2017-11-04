
from numpy import arange
import math
from matplotlib.pyplot import *

g=6.673*10**-8
c=299792458
e0= 1
R0 = 1.47
K = 5
powerstuff = 4/3.
p0 = 1000

alpha = 1.476
beta = 0.03788
#alpha = R0/(K*(e0**(powerstuff-1)))**(1/powerstuff)
#beta = (4*math.pi*e0)/((c**2)*(K*e0**(powerstuff-1))**(1/powerstuff))
# BOTH OF THE ABOVE ARE IN M^2, CHANGE TO KM AT SOME POINT
# need some function for density

#NEWTONIAN NON-RELATIVISTIC WHITE DWARF
# SOLAR MASSES & UNITLESS
#
def diffpressure(p, r, m):
     return -(alpha*(p**(1./powerstuff))*m/(r**2))
#    return -(g*epsilon(r)*mass(r))/((r*c)**2)
def diffmass(p, r):
    return (p**(1./powerstuff))*beta*r**2




def solve_rk4_coupled(mass, pressure, p0, m0, N, rinitial, rfinal):
    p = p0
    m = m0
    h = (rfinal - rinitial) / float(N)
    rs = arange(rinitial, rfinal, h)
    ps = []
    ms = []
    for r in rs:
        print r, m, p
        ms.append(m)
        ps.append(p)
        k1 = h * mass(p, r)
        l1 = h * pressure(p, r, m)
        k2 = h * mass(p + k1 / 2., r + h / 2.)
        l2 = h * pressure(p + l1 / 2., r + h / 2., m + k1 / 2.)
        k3 = h * mass(p + k2 / 2., r + h / 2.)
        l3 = h * pressure(p + l2 / 2., r + h / 2., m + k2 / 2.)
        k4 = h * mass(p + k3, r + h)
        l4 = h * pressure(p + l3, r + h, m + k3)
        p += 1 / 6. * (k1 + 2 * k2 + 2 * k3 + k4)
        m += 1 / 6. * (l1 + 2 * l2 + 2 * l3 + l4)
    return rs, ps, ms

rs, ps, ms = solve_rk4_coupled(diffmass, diffpressure, p0, 100, 100, 0.1, 10000)
plot(rs, ps, 'r')
#plot(rs, ms, 'b')


show()