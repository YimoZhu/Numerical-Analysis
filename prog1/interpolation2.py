# -*- coding: utf-8 -*-
"""
Created on Fri Oct 05 16:42:21 2018

@author: 49048
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from interpolation1 import output1

def barycentricInterpolation(Xs,fXs,e):
    pairs = dict(zip(Xs,fXs))
    if e in Xs:
        result = pairs[e]
    else:
        j = 0
        numerator = 0
        denominator = 0
        for xj in Xs:
            #calculate lambda j
            k = 0
            lambda_j_inv = 1
            for xk in Xs:
                if k != j:
                    lambda_j_inv = lambda_j_inv*(xj-xk)
                else:
                    pass
                k = k+1
            lambda_j = 1/lambda_j_inv
            numerator = numerator + lambda_j*fXs[j]/(e-xj)
            denominator = denominator + lambda_j/(e-xj)
            j = j+1
        result = numerator/denominator
    return result

def output2(f,xin,xout):
    Xs = xin
    fXs = [f(x) for x in Xs]
    yout = np.array([barycentricInterpolation(list(Xs),fXs,x) for x in xout])
    return yout

def main():
    Ns = [10,19,50,99]
    f = lambda x:1.0/(1+9*x**2)
    xout = np.linspace(-1,1,500)
    for n in Ns:
        xin1 = np.linspace(-1,1,n+1)
        xin2 = np.cos(np.linspace(-np.pi,0,n+1))
        print '______________n = %s_________________'%n
        print output2(f,xin1,xout)
        print output2(f,xin2,xout)
        print ''
        
main()

"""The following codes are to generate the 5-8 graphics"""
f = lambda x:1.0/(1+9*x**2)
xout = np.linspace(-1,1,500)
y = np.array([f(x) for x in xout])
#Case 1: n = 50, xin1
n = 50
xin1 = np.linspace(-1,1,n+1)
yout1 = output1(f,xin1,xout)
yout2 = output2(f,xin1,xout)
fig5 = plt.figure(figsize=(10,8))
fig5.suptitle('Semilogy plot, n = 50 , uniform interpolation')
plt.semilogy(xout,1.0e-18+abs(yout1-y),label='Lagrange')
plt.semilogy(xout,1.0e-18+abs(yout2-y),label='barycentric')
plt.legend()
plt.xlim=[-1,1]
plt.show()
#Case 2: n = 50, xin2
n = 50
xin1 = np.linspace(-1,1,n+1)
xin2 = np.cos(np.linspace(-np.pi,0,n+1))
yout1 = output1(f,xin2,xout)
yout2 = output2(f,xin2,xout)
fig6 = plt.figure(figsize=(10,8))
fig6.suptitle('Semilogy plot, n = 50 , Chebyshev interpolation')
plt.semilogy(xout,1.0e-18+abs(yout1-y),label='Lagrange')
plt.semilogy(xout,1.0e-18+abs(yout2-y),label='barycentric')
plt.legend()
plt.show()
#Case 3: n = 99, xin1
n = 99
xin1 = np.linspace(-1,1,n+1)
xin2 = np.cos(np.linspace(-np.pi,0,n+1))
yout1 = output1(f,xin1,xout)
yout2 = output2(f,xin1,xout)
fig7 = plt.figure(figsize=(10,8))
fig7.suptitle('Semilogy plot, n = 99 , Uniform interpolation')
plt.semilogy(xout,1.0e-18+abs(yout1-y),label='Lagrange')
plt.semilogy(xout,1.0e-18+abs(yout2-y),label='barycentric')
plt.legend()
plt.show()
#Case 4: n = 99, xin2
n = 99
xin1 = np.linspace(-1,1,n+1)
xin2 = np.cos(np.linspace(-np.pi,0,n+1))
yout1 = output1(f,xin2,xout)
yout2 = output2(f,xin2,xout)
fig8 = plt.figure(figsize=(10,8))
fig8.suptitle('Semilogy plot, n = 99 , Chebyshev interpolation')
plt.semilogy(xout,1.0e-18+abs(yout1-y),label='Lagrange')
plt.semilogy(xout,1.0e-18+abs(yout2-y),label='barycentric')
plt.legend()
plt.show()