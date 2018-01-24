# NEWTONIAN RELATIVISTIC WHITE DWARF
# SOLAR MASSES & UNITLESS
# Kaitlin Williams
#
#
#
import numpy as np
import math
from matplotlib.pyplot import *

#seterr(all='print')

g=6.673*(10**-8)
c=299792458
solarmass = 1.989*(10**30)
powerstuff = 5./3
R0 = 1.473
K = (1./4.173)**(powerstuff-1)

p0 = [10**-4]
#p0 = np.arange(0.001,1.00,0.001,object)


#-------------------------------
# EQUATION OF STATE
#alpha = R0/((K*(e0**(powerstuff-1))))**(1/powerstuff)


e0= 0.08969
alpha = 1.
beta = 0.7636

#beta = (4*math.pi*e0)/((c**2)*solarmass*(K*(e0**(powerstuff-1)))**(1/powerstuff))

#-------------------------------



def diffpressure(p, r, m):
    '''
    Diffpressure acts as the Newtonian differential pressure equation (one of two coupled equations involved in the TOV).
    If pressure is negative, throws large number wrench in the works to stop the process.

    :param p: pressure
    :param r: radius
    :param m: mass
    :return: numerical output of equation 23
    '''

     if m ==0:
        return p
     if p <= 0 :
        print "I'M NEGATIVE  ", p
        return -1
     else:
        return -(alpha*m*(p**(1./powerstuff)))/(r**2)


def diffmass(p, r):
    '''
    Diffmass acts as the Newtonian differential mass equation (the second coupled equation involved in the TOV).
    If pressure is negative, throws large number wrench in the works to stop the process.

    :param p: pressure
    :param r: radius
    :return: numerical output of equation 26
    '''
    if p <= 0 :
        print "I'M NEGATIVE ", p
        return -1
    return (p**(1./powerstuff))*beta*(r**2)




def solve_rk4_coupled(mass, pressure, p0, m0, N, rinitial, rfinal):
    '''
    Coupled RK4 equation solver.
    :param mass: differential mass equation
    :param pressure: differential pressure equation
    :param p0: initial pressure
    :param m0: initial mass (currently always 0 but this could change)
    :param N: number of iterations wanted
    :param rinitial: initial radius (currently always 0 but this could change)
    :param rfinal: final radius reached by N iterations. Usually around 20-30.
    :return: returns either a list of 1x3 arrays, or NULL if the program cannot find when p reaches a negative value in
    the r range of rinitial-rfinal.

    '''
    p = p0
    m = m0
    r = rinitial
    h = (rfinal - rinitial) / float(N)
    data = [[r, p, m]]
    for i in range(1,N):
        if i % 10000 == 0:
            print ".",
        l1 = h * mass(p, r)
        k1 = h * pressure(p, r, m)
        l2 = h * mass(p + k1 / 2., r + h / 2.)
        k2 = h * pressure(p + k1 / 2., r + h / 2., m + l1 / 2.)
        l3 = h * mass(p + k2 / 2., r + h / 2.)
        k3 = h * pressure(p + k2 / 2., r + h / 2., m + l2 / 2.)
        l4 = h * mass(p + k3, r + h)
        k4 = h * pressure(p + k3, r + h, m + l3)
        r += h
        p += (k1 + (2 * k2) + (2 * k3) + k4)/6.
        m += (l1 + (2 * l2) + (2 * l3) + l4)/6.
        if p <= 0:
            break
        data.append([r, p, m])
    else:
        print data[-1]
        data = None

    return data

finalrs = []
finalms = []


for press in p0:
    stardata = solve_rk4_coupled(diffmass, diffpressure, press, 0., 10000001, 0, 30.)
    rs = []
    ms = []
    ps = []

    for elem in stardata:
        rs.append(elem[0])
        ps.append(elem[1])
        ms.append(elem[2])
#    finalrs.append(stardata[-1][0])
#    finalms.append(stardata[-1][2])
#    print stardata[-1]

    plot(rs, ps, 'r')
    #plot(rs, ms, 'b')

#plot(finalrs, finalms, 'r')
show()
