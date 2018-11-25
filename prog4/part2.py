from numpy import *
from rkf import rkf
from matplotlib import pyplot as plt
from sdirk import sdirk

def f(y):
    #The function handel
    return array([-.04*y[0]+y[1]*y[2],
                  400*y[0]-1e4*y[1]*y[2]-3000*y[1]**2,
                  .3*y[1]**2])

def J(y):
    #The Jacobian matrix of f evaluated at point y.
    return array([[-.04,y[2],y[1]],
                  [400,-1e4*y[2]-6000*y[1],-1e4*y[1]],
                  [0,.6*y[1],0]])

t2,y2 = rkf(f,0,3,[1,0,0],0.1,1e-7,0.25,1e-3)
t3,y3 = sdirk(f,J,0,3,[1,0,0],0.1,1e-7,0.25,1e-3)
fig = plt.figure(figsize=(6,4))
plt.suptitle("RKF method")
plt.plot(t2,y2.T)
plt.grid(True)
fig = plt.figure(figsize=(6,4))
plt.suptitle("SDIRK method")
plt.plot(t3,y3.T)
plt.grid(True)
plt.show()