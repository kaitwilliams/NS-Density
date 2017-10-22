import math
a = input("Please give the value of a ")
b = input("Please give the value of b ")
c = input("PLease give the value of c ")
bottom = 2*a
first_top = -b + (b**2 - 4*a*c)**0.5
second_top = -b - (b**2 - 4*a*c)**0.5
x1 = float(first_top)/bottom
x2 = float(second_top)/bottom
print "The values of x are " + str(x1) + str(x2)
deriv_first_bottom = b + (b**2 - 4*a*c)**0.5
deriv_second_bottom = b - (b**2 - 4*a*c)**0.5
deriv_top=-2*c
d_x1 = float(deriv_top)/deriv_first_bottom
d_x2 = float(deriv_top)/deriv_second_bottom
print "The values of x' are " + str(d_x1) + str(d_x2)
