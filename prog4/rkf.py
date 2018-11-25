from numpy import *

def rkf(f,a,b,y0,h,hmin,hmax,tol):
    t=[a]
    y = [y0]
    h = min(h,b-a)

    A=array([[0,0,0,0,0,0],
             [1/4,0,0,0,0,0],
             [3/32,9/32,0,0,0,0],
             [1932/2197,-7200/2197,7296/2197,0,0,0],
             [439/216,-8,3680/513,-845/4104,0,0],
             [-8/27,2,-3544/2565,1859/4104,-11/40,0]])

    bb = array([25/216,          0,  1408/2565,   2197/4104,   -1/5,    0 ])
    b1 = array([  16/135,          0, 6656/12825, 28561/56430,  -9/50, 2/55 ])
    er = bb - b1
    
    s = len(b1)
    d = len(y0)
    kk = zeros((d,s))
    
    nmax = 10000
    for n in range(1,nmax+1):
        for i in arange(s):
            z = array(y0)
            for j in arange(i):
                #Compute the stage derivatives
                z = z+h*A[i,j]*kk[:,j]
            kk[:,i] = f(z)

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
    
    print("[RKF] Iteration ended at loop # %s"%n)
    return (array(t),array(y).T)
        