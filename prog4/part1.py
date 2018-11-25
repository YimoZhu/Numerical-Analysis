from newton import newton
from numpy import *

#Testing the newton solver.
def f(x):
    #x is a 4 dim vector.
    return x-exp(cos(arange(1,5)*sum(x)))

def J(x):
    #x is a 4 dim vector, return f's Jacobian matrix evaluated at x.
    J = eye(4)
    u = sum(x)
    for i in range(4):
        for j in range(4):
            J[i,j] = J[i,j] - (i+1)*exp((i+1)*u)
    return J

#Then solve the equation f(x) = 0
p1 = newton(f,J,array([2.5,2,1.4,0.9]),1e-12)
print(p1)
print(f(p1[0]))