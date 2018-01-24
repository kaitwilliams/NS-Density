# NEWTONIAN NON-RELATIVISTIC WHITE DWARF
# SOLAR MASSES & UNITLESS
# Kaitlin Williams
#
#
#
import numpy as np
import math
from matplotlib.pyplot import *

#seterr(all='print')

g=6.673*10**-8
c=299792458
#solarmass =
powerstuff = 4/3.
e0= 4.173
R0 = 1.473
K = (1./4.173)**(powerstuff-1)
p0 = [1.00]
#p0 = np.arange(0.001,1.00,0.001,object)


#-------------------------------

alpha = R0/((K*(e0**(powerstuff-1))))**(1/powerstuff)
#print alpha
#alpha = 1.473
beta = 52.46

#beta = (4*math.pi*e0)/((c**2)*(K*(e0**(powerstuff-1)))**(1/powerstuff))
# ^^ needs the solar mass term otherwise waaaay too small
#print beta
#-------------------------------


def diffpressure(p, r, m):
     if m ==0:
        return p
     if p <=0:
 #       print "I'M NEGATIVE"
        return -1
     else:
        return -(alpha*m*(p**(1./powerstuff)))/(r**2)
#    return -(g*epsilon(r)*mass(r))/((r*c)**2)

def diffmass(p, r):
    return (p**(1./powerstuff))*beta*r**2




def solve_rk4_coupled(mass, pressure, p0, m0, N, rinitial, rfinal):
    '''

    :param mass:
    :param pressure:
    :param p0:
    :param m0:
    :param N: number of iterations
    :param rinitial: starting radius
    :param rfinal:
    :return:
    '''
    p = p0
    m = m0
    r = rinitial
    h = (rfinal - rinitial) / float(N)
    data = [[r, p, m]]
    for i in range(1,N):
        r = r+h
        k1 = h * mass(p, r)
        l1 = h * pressure(p, r, m)
        k2 = h * mass(p + k1 / 2., r + h / 2.)
        l2 = h * pressure(p + l1 / 2., r + h / 2., m + k1 / 2.)
        k3 = h * mass(p + k2 / 2., r + h / 2.)
        l3 = h * pressure(p + l2 / 2., r + h / 2., m + k2 / 2.)
        k4 = h * mass(p + k3, r + h)
        l4 = h * pressure(p + l3, r + h, m + k3)
        p += 1 / 6. * (l1 + 2 * l2 + 2 * l3 + l4)
        m += 1 / 6. * (k1 + 2 * k2 + 2 * k3 + k4)
        if p <= 0:
            break
        data.append([r, p, m])
    else:
 #       print data[-1]
        data = None

    return data

finalrs = []
finalms = []


for press in p0:
    stardata = solve_rk4_coupled(diffmass, diffpressure, press, 0., 1000000, 0, 1000.)
    rs = []
    ms = []
    ps = []

    for elem in stardata:
        rs.append(elem[0])
        ps.append(elem[1])
        ms.append(elem[2])
    finalrs.append(stardata[-1][0])
    finalms.append(stardata[-1][2])
    print stardata[-1]
    plot(rs, ps, 'r')
    #plot(rs, ms, 'b')

#(finalrs, finalms, 'r')
show()
