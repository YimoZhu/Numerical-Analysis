# -*- coding: utf-8 -*-
"""
Created on Fri Oct 05 16:20:59 2018

@author: 49048
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt

def LagrangeInterpolate(Xs,fXs,e):
    result = 0
    j = 0
    for xj in Xs:
        Ljx = 1
        k = 0
        for xk in Xs:
            if k != j:  
                Ljx = Ljx*(e-xk)/(xj-xk)
            else:
                pass
            k = k+1
        result = result + fXs[j]*Ljx
        j = j+1
    return result

def output1(f,xin,xout):
    Xs = xin
    fXs = [f(x) for x in Xs]
    yout = np.array([LagrangeInterpolate(Xs,fXs,x) for x in xout])
    return yout

def main():
    Ns = [10,19,50,99]
    f = lambda x:1.0/(1+9*x**2)
    xout = np.linspace(-1,1,500)
    for n in Ns:
        xin1 = np.linspace(-1,1,n+1)
        xin2 = np.cos(np.linspace(-np.pi,0,n+1))
        print '______________n = %s_________________'%n
        print output1(f,xin1,xout)
        print output1(f,xin2,xout)
        print ''

main()

"""The following codes are to generate the 1-4 graphics"""
#plot fx & px when n = 10
f = lambda x:1.0/(1+9*x**2)
xout = np.linspace(-1,1,500)
y = np.array([f(x) for x in xout])
fig1 = plt.figure(figsize=(10,8))
fig1.suptitle('n=10 & using Lagrange P(x) & uniform')
plt.plot(xout,y,label='f(x)')
n = 10
xin1 = np.linspace(-1,1,n+1)
xin2 = np.cos(np.linspace(-np.pi,0,n+1))
px = output1(f,xin1,xout)
plt.plot(xout,px,label=r'$P_{10}(x)$')
plt.legend()
plt.show()
#plot fx & px when n = 10
fig2 = plt.figure(figsize=(10,8))
fig2.suptitle('n = 10 & using Lagrange P(x) & Chebyshev')
px = output1(f,xin2,xout)
plt.plot(xout,y,label='f(x)')
plt.plot(xout,px,label=r'$P_{10}(x)$')
plt.legend()
plt.show()
#plot fx & px when n = 19
fig3 = plt.figure(figsize=(10,8))
fig3.suptitle('n=19 & using Lagrange P(x) & uniform')
plt.plot(xout,y,label='f(x)')
n = 19
xin1 = np.linspace(-1,1,n+1)
xin2 = np.cos(np.linspace(-np.pi,0,n+1))
px = output1(f,xin1,xout)
plt.plot(xout,px,label=r'$P_{19}(x)$')
plt.legend()
plt.show()
#plot fx & px when n = 19
fig4 = plt.figure(figsize=(10,8))
fig4.suptitle('n = 19 & using Lagrange P(x) & Chebyshev')
px = output1(f,xin2,xout)
plt.plot(xout,y,label='f(x)')
plt.plot(xout,px,label=r'$P_{19}(x)$')
plt.legend()
plt.show()
