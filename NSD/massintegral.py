from numpy import arange
import math
from matplotlib.pyplot import *


def f(x, t):
    return -x ** 4 + math.sin(2. * t)


def solve_euler(f, a, b, N, xinit):
    x = xinit
    h = (b - a) / float(N)
    ts = arange(a, b, h)
    xs = []
    for t in ts:
        xs.append(x)
        x += h * f(x, t)
    return ts, xs


def solve_rk2(f, a, b, N, xinit):
    x = xinit
    h = (b - a) / float(N)
    ts = arange(a, b, h)
    xs = []
    for t in ts:
        xs.append(x)
        k1 = h * f(x, t)
        k2 = h * f(x + k1 / 2., t + h / 2.)
        x += k2
    return ts, xs


def solve_rk4(f, a, b, N, xinit):
    x = xinit
    h = (b - a) / float(N)
    ts = arange(a, b, h)
    xs = []
    for t in ts:
        xs.append(x)
        k1 = h * f(x, t)
        k2 = h * f(x + k1 / 2., t + h / 2.)
        k3 = h * f(x + k2 / 2., t + h / 2.)
        k4 = h * f(x + k3, t + h)
        x += 1 / 6. * (k1 + 2 * k2 + 2 * k3 + k4)
    return ts, xs


ts, xs = solve_rk4(f, 0., 10., 20, 0.0)
plot(ts, xs, 'r')
ts, xs = solve_rk4(f, 0., 10., 50, 0.0)
plot(ts, xs, 'g')
ts, xs = solve_rk4(f, 0., 10., 100, 0.0)
plot(ts, xs, 'b')
ts, xs = solve_rk4(f, 0., 10., 1000, 0.0)
plot(ts, xs, 'y')

show()

