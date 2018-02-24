
# NEWTONIAN RELATIVISTIC WHITE DWARF
# SOLAR MASSES & UNITLESS
# Kaitlin Williams
#
#
#

import numpy as np
import math
from matplotlib.pyplot import *

np.seterr(all='print')

g=6.673*(10**-8)
c=299792458
solarmass = 1.989*(10**30)

R0 = 1.476

p0 = [0.0562341]
i=0.00005

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


e0= 0.08969
alpha = 1.476
beta = 0.03778

#beta = (4*math.pi*e0)/((c**2)*solarmass*(K*(e0**(powerstuff-1)))**(1/powerstuff))

#-------------------------------

def eos(p):
    return (aNR*(p**powerstuff1) + aR*(p**powerstuff2))


def diffpressureNewtonian(p,r,m):
    '''
     Diffpressure acts as the Newtonian differential pressure equation (one of two coupled equations involved in the TOV).
     If pressure is negative, throws large number wrench in the works to stop the process.
     :param p: pressure
     :param r: radius
     :param m: mass
     :return: numerical output of equation 23
     '''

    if m == 0:
        return p
    if p <= 0:
        return -1
    else:
        #new_EOS = (aNR * (p ** powerstuff1) + aR * (p ** powerstuff2))
        term0 = -(alpha * m * (eos(p))) / (r ** 2)
        return term0


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
        return -1
    else:
       # new_EOS = (aNR*(p**powerstuff1) + aR*(p**powerstuff2))
        term0 = -(alpha * m * (eos(p))) / (r ** 2)
        term1 = (1 + (p/(eos(p))) * (R0 / alpha))
        term2 = (1 + ((beta * R0) / (m * alpha)) * p * (r ** 3))
        term3 = (1 - (2 * R0 * m) / r) ** (-1)
        return term0*term1*term2*term3


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
    #new_EOS = (aNR * (p ** powerstuff1) + aR * (p ** powerstuff2))

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
      #  if i % 100000 == 0:
          #  print ".",
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
        data = None

    return data

finalrs = []
finalms = []


for press in p0:
    #stardata = solve_rk4_coupled(diffmass, diffpressure, press, 0., 1000, 0, 20)
    stardata2 = solve_rk4_coupled(diffmass, diffpressureNewtonian, press, 0., 1000, 0, 20)
    rs = []
    rs2 = []
    ms = []
    ms2 = []
    ps = []
    ps2 = []

    for elem in stardata2:
        rs.append(elem[0])
        ps.append(elem[1])
        ms.append(elem[2])
  #  finalrs.append(stardata2[-1][0])
   # finalms.append(stardata2[-1][2])
    print stardata2[-1]

    plot(rs, ms, 'b')

  #  for elemn in stardata2:
   #     rs2.append(elemn[0])
   #     ps2.append(elemn[1])
   #     ms2.append(elemn[2])

    #plot(rs2, ms2, 'r')

    #plot(rs, ms, 'b')

#plot(finalrs, finalms, 'r')
show()