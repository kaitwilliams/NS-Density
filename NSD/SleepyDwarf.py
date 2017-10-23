
from numpy import arange
import math
from matplotlib.pyplot import *

#need definitions here so can use in equations
#g=6.673*10^-8
#c=speed of light
#need some function for density

#def f(p,r):
#	return
#def g(m,r):
#	return 4pir^2epsilon(r)/c^2 <--ooo that'll be trouble

#unedited RK2 RK4 code. Need to modify for range of numbers and 
incorporate second function.

def solve_rk2(f,a,b,N,xinit):
	x = xinit
	h = (b-a)/float(N)
	ts = arange(a,b,h)
	xs = []
	for t in ts:
		xs.append(x)
		k1 = h*f(x,t)
		k2 = h*f(x+k1/2., t+h/2.)
		x += k2
	return ts,xs

def solve_rk4(f,a,b,N,xinit):
	x = xinit
	h = (b-a)/float(N)
	ts = arange(a,b,h)
	xs = []
	for t in ts:
		xs.append(x)
		k1 = h*f(x,t)
		k2 = h*f(x+k1/2., t+h/2.)
		k3 = h*f(x+k2/2., t+h/2.)
		k4 = h*f(x+k3,t+h)
		x += 1/6.*(k1+2*k2+2*k3+k4)
	return ts,xs

 
