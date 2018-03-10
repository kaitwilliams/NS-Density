# NEWTONIAN RELATIVISTIC WHITE DWARF
# SOLAR MASSES & UNITLESS
# Kaitlin Williams
#
#
#

import numpy as np
import math
from scipy.integrate import odeint
from matplotlib.pyplot import *

np.seterr(all='print')

g=6.673*(10**-8)
c=2.998*(10**10)
solarmass = 1.989*(10**33)
hbar = 1.055*(10**-27)
mneutron = 1.67*(10**-24)

km3 = 10**15
Mc2 = solarmass*(c**2)
#e0init = (mneutron**4)*(c**8)/(3*(math.pi**2)*(hbar**3)*(c**3))
#e0 = e0init*(km3/Mc2)
e0 = 0.003006

R0 = (g * solarmass/(c**2))/10**5
#p0 = [0.01]
p0 = np.logspace(-4,3,20)

#while (i<20.0):
 #   p0 += [i]
  #  i += 0.001
#p0 = [0.01]

#-------------------------------
# EQUATION OF STATE
# alpha = R0/((K*(e0**(powerstuff-1))))**(1/powerstuff)

powerstuff1 = 3./5
aNR = 2.4216
powerstuff2 = 1.
aR = 2.8663
#K = (1./4.173)**(powerstuff-1)



#alpha = R0
alpha = 1.476
#beta = (4*3.14*e0)
beta = 0.03778
#beta = (4*math.pi*e0)/((c**2)*solarmass*(K*(e0**(powerstuff-1)))**(1/powerstuff))

#-------------------------------

def eos(p):
    return (aNR*(p**powerstuff1) + aR*(p**powerstuff2))


def diffpressure(p, r, m, p0):
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
        return -1
    else:
        #new_EOS = (aNR*(p**powerstuff1) + aR*(p**powerstuff2))
#        term0 = -(alpha * m * (eos(p))) / (r ** 2)
        term0 = - (alpha * eos(p) * m)/(r**2)
#        term1 = (1 + (p/(eos(p)))*(R0 / alpha))
        term1 = 1 + (p/eos(p))
   #     term2 = (1 + ((beta * R0 * p * (r ** 3))/ (m * alpha)))

        if r < 0.1:
            term2 = 1 + ((3*e0)/(solarmass * c**2))*p/eos(p0)
        else:
            term2 = 1 + ((4 * math.pi * e0)/(solarmass * c**2)) * r**3 * p/m
  #      term3 = (1 - (2 * R0 * m)/r)**(-1)

        term3 = 1 - (2*R0*m)/r
        #term3 = 1 - (2*g*solarmass*m)/((c**2)*r)
        return (term0*term1*term2)/term3


def diffmass(p, r):
    '''
    Diffmass acts as the Newtonian differential mass equation (the second coupled equation involved in the TOV).
    If pressure is negative, throws large number wrench in the works to stop the process.
    :param p: pressure
    :param r: radius
    :return: numerical output of equation 26
    '''
    if p <= 0 :
        return -1
   # new_EOS = (aNR * (p ** powerstuff1) + aR * (p ** powerstuff2))

    return (eos(p))*beta*(r**2)




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
        l1 = h * mass(p, r)
        k1 = h * pressure(p, r, m, p0)

        l2 = h * mass(p + k1 / 2., r + h / 2.)
        k2 = h * pressure(p + k1 / 2., r + h / 2., m + l1 / 2., p0)

        l3 = h * mass(p + k2 / 2., r + h / 2.)
        k3 = h * pressure(p + k2 / 2., r + h / 2., m + l2 / 2., p0)

        l4 = h * mass(p + k3, r + h)
        k4 = h * pressure(p + k3, r + h, m + l3, p0)

        r += h
        p += (k1 + (2 * k2) + (2 * k3) + k4)/6.
        m += (l1 + (2 * l2) + (2 * l3) + l4)/6.
        if p <= 0:
            break
        data.append([r, p, m])
    else:
        data = None

    return data

finalrs = []
finalms = []


for press in p0:
    stardata = solve_rk4_coupled(diffmass, diffpressure, press, 0., 2000, 0, 30)
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


    #plot(rs, ps, 'r')
    #plot(rs, ms, 'b')

plot(finalrs, finalms, 'r')
show()

#print stardata[-1]