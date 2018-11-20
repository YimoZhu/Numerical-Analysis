from numpy import *
from rkf import rkf
from matplotlib import pyplot as plt

def f(y):
    #The function handel
    return array([-.04*y[0]+y[1]*y[2],
                  400*y[0]-1e4*y[1]*y[2]-3000*y[1]**2,
                  .3*y[1]**2])

def J(y):
    #The Jacobian matrix of f evaluated at point y.

t2,y2 = rkf(f,0,3,[1,0,0],0.1,1e-7,0.25,1e-3)
fig = plt.figure(figsize=(6,4))
plt.plot(t2,y2)
