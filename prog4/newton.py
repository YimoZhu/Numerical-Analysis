from numpy import *

def newton(f,J,p0,tol):
    #input the function and its Jacobian matrix J, and the starting point p0 to get the zero of f.
    f0 = f(p0)
    J0 = J(p0)
    pp = array([[1,linalg.norm(f0),linalg.norm(J0)]])
    for n in arange(50):
        p = p0 - linalg.inv(J(p0)).dot(f(p0))
        f0 = f(p)
        J0 = J(p)
        pp = append(pp,array([[n+2,linalg.norm(f0),linalg.norm(J0)]]),axis=0)
        if linalg.norm(p-p0) < tol:
            break
        else:
            p0 = p
            continue
    return (p,n,pp)