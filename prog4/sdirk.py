from numpy import *
from newton import newton

def sdirk(f,J,a,b,y0,h,hmin,hmax,tol):
    t=[a]
    y = [y0]
    h = min(h,b-a)

    A=array([[1/4,0,0,0,0],
             [1/2,1/4,0,0,0],
             [17/50,-1/25,1/4,0,0],
             [371/1360,-137/2720,15/544,1/4,0],
             [25/24,-49/48,125/16,-85/12,1/4]])

    bb = array([25/24,          -49/48,  125/16,   -85/12,   1/4])
    b1 = array([59/48,-17/96,225/32,-85/12,0])
    er = bb - b1
    
    s = len(b1)
    d = len(y0)
    kk = zeros((d,s))
    
    nmax = 10000
    for n in range(1,nmax+1):
        for i in arange(s):
            #Calculate the stage derivative.
            z = array(y0)
            for j in arange(i):
                z = z+h*A[i,j]*kk[:,j]
            ah = h*A[i,i]
            g = lambda k:k-f(z+ah*k)
            dg = lambda k:eye(d)-J(z+ah*k)*ah
            kk[:,i] = newton(g,dg,f(y0),1.0e-12)[0]

        R = linalg.norm(kk.dot(er))

        if R<tol:
            #Accept step
            t.append(t[-1]+h)
            w = list(array(y0) + h*kk.dot(bb))
            y.append(w)
            y0 = w
            if t[-1] > a+0.999999999*(b-a):
                break

        q = (0.5*tol/R)**0.25
        q = max(q,0.1)
        q = min(q,4.0)
        h = q*h
        if h <hmin:
            h = hmin
        if h >hmax:
            h = hmax
        if h > b-t[-1]:
            h = b - t[-1]
    if n >= nmax:
        print("error: nmax reached %s"%s)
    
    print("[SDIRK] Iteration ended at loop # %s"%n)
    return (array(t),array(y).T)
        