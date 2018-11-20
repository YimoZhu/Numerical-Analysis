from numpy import *

def newton(f,J,p0,tol):
    #input the function and its Jacobian matrix J, and the starting point p0 to get the zero of f.
    for i in arange(1,101):
        p = p0 - linalg.inv(J(p0)).dot(f(p0))
        if linalg.norm(p-po) < tol:
            break
        else:
            p0 = p
            continue
    return p