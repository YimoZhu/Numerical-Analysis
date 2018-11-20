# -*- coding: utf-8 -*-
"""
@Name:some numerical algorithms.
@author: Yimo Zhu
"""
from __future__ import division
import numpy as np
from matplotlib import pyplot as plt

def aitken(g,p0,steps=20,TOL=1e-5,mode='steps'):
    if mode == 'steps':
        sequence = []
        sequence.append(p0)
        phat = []
        print 'calculating original sequence...'
        print 'p0 = %s'%p0
        for i in range(steps+1):
            p = g(p0)
            sequence.append(p)
            print 'p%s = %s'%(i+1,p)
            p0 = p
        print ''
        print 'calculating accelerated sequence...'
        for i in range(steps):
            pn = sequence[i]
            pn1 = sequence[i+1]
            pn2 = sequence[i+2]
            delta = pn1-pn
            delta2 = pn2-2*pn1+pn
            result = pn-delta**2/delta2
            print 'phat%s = %s'%(i,result)
            phat.append(result)
        return type('result_aitken',(),{'original_seq':sequence,'phat':phat})
    
def steffensen(g,p0,steps=100000,TOL=0,mode='steps'):
    phats = []
    phat = None
    for i in range(steps):
        p1 = g(p0)
        p2 = g(p1)
        delta = p1 - p0
        delta2 = p2 - 2*p1 + p0
        phat_new = p0 - delta**2/delta2
        phats.append(phat_new)
        print 'phat%s = %s'%(i,phat_new)
        if phat == None:
            pass
        else:
            dev = abs(phat_new-phat)
            if dev < TOL:
                break
            else:
                pass
        p0 = phat_new
        phat = phat_new
    return type('result_steffensen',(),{'phat':phats})

def newton(f,f_prime,p0,steps=10000,TOL=1e-5):
    g = lambda x:(x-f(x)/f_prime(x))
    sequence = []
    sequence.append(p0)
    print 'p0 = %s'%p0
    for i in range(steps-1):
        p = g(p0)
        sequence.append(p)
        print 'p%s = %s'%(i+1,p)
        if abs(p-p0)<TOL:
            break
        p0 = p
    return type('result_newtown',(),{'p':sequence})

def plot(f,interval):
    x = np.linspace(interval[0],interval[1],256)
    y = np.array([f(i) for i in x])
    plt.figure()
    plt.plot(x,y)
    plt.grid()
    
def homer(coef,x0):
    coef1 = []
    coef1.append(coef[0])
    for i in range(len(coef)-1):
        coef1.append(coef[i+1]+x0*coef1[-1])
    b0 = coef1.pop()
    coef2 = []
    coef2.append(coef1[0])
    for i in range(len(coef1)-1):
        coef2.append(coef1[i+1]+x0*coef2[-1])
    prime=coef2.pop()
    return type('result_homer',(),{'b0':b0,'derivative':prime,'Qx':coef1,'residual':coef2})

class polynomial_deg3(object):
    def __init__(self):
        #Constructor,no other usage,
        #just to clarify the important parameters of a\
        #polynomial with degree up to 3
        self.c=None
        self.d=None
        self.e=None
        self.f=None
        self.degree=None
        
    def set_parameter(self,iterable):
        #To set the parameter c,d,e,f of the equation,using a iterable object        
        self.c,self.d,self.e,self.f = tuple(iterable)
        #To determine the degree of the polynomial,
        #according to the parameters,in convenience \
        #for solving the root of derivative
        degree=4
        for coef in iterable:
            degree = degree-1
            if coef != 0:
                break
            else:
                continue
        self.degree = degree
    
    def evaluate(self,x):
        #To evaluate the polynomial at x
        return self.c*x**3+self.d*x**2+self.e*x+self.f

    def evaluate_derivative(self,x):
        #To evaluate the derivative of the polynomial at x
        return 3*self.c*x**2 + 2*self.d*x+self.e
    
    def get_derivative_solution(self):
        #To solve p'(x)=0
        if self.degree==3:
            #If the polynomial has degree3, then p'(x)=0 is a quadratic equation
            discriminant = (2*self.d)**2 - 4*(3*self.c)*(self.e)
            if discriminant < 0:
                return 'No root!'
            elif discriminant == 0:
                return -(2*self.d)/(2*(3*self.c))
            else:
                #This is using proper quadratic formula
                #case 1: if d and the discriminant has different signal
                if self.d > 0:
                    x1 = (-2*self.d-sqrt(discriminant))/(2*3*self.c)
                    x2 = (2*self.e)/(-2*self.d-sqrt(discriminant))
                    return (x1,x2)
                elif self.d < 0:
                    x1 = (-2*self.d+sqrt(discriminant))/(2*3*self.c)
                    x2 = (2*self.e)/(-2*self.d+sqrt(discriminant))
                    return (x1,x2)
                else:#self.d == 0
                    return (sqrt(-self.e/(3*self.c)),-sqrt(-self.e/(3*self.c)))
        elif self.degree==2:
            #If the polynomial has degree2, then p'(x)=0 is a linear eaquation
            return -(self.e)/(2*self.d)
        elif self.degree==1:
            return 'No root!'
        else:
            return 'No root!'

def muller(f,p0,p1,p2,steps=10000,TOL=1e-4):
    sequence = [p0,p1,p2]
    for i in range(steps-3):
        a = ((p1-p2)*(f(p0)-f(p2))-(p0-p2)*(f(p1)-f(p2)))/((p0-p2)*(p1-p2)*(p0-p1))
        b = ((p0-p2)**2*(f(p1)-f(p2))-(p1-p2)**2*(f(p0)-f(p2)))/((p0-p2)*(p1-p2)*(p0-p1))
        c = f(p2)
        det = b**2-4*a*c
        if type(det) == complex:
            rootdet = np.sqrt(det)
        elif det >= 0:
            rootdet=np.sqrt(det)
        elif det<0 :
            rootdet = np.sqrt(-det)*1j
        p3 = p2 - (2*c)/(b+np.sign(b)*rootdet)
        sequence.append(p3)
        print 'p%s = %s'%(i+3,p3)
        if abs(p3-p2)<TOL:
            break
        p0 = p1
        p1 = p2
        p2 = p3
    return sequence


"""
Interpolating related methods
"""
def neville(Xs,fXs,evaluation_point):
    if len(Xs) != len(fXs):
        raise ValueError
    n = len(Xs) - 1
    chart = np.array([[np.nan]*(n+1)]*(n+1))
    chart[:,0] = fXs
    for i in 1+np.arange(n):
        for j in 1+np.arange(i):
            chart[i,j] = ((evaluation_point-Xs[i-j])*chart[i,j-1]-(evaluation_point-Xs[i])*chart[i-1,j-1])/(Xs[i] - Xs[i-j])
    print chart
    return chart

def compositeSimpson(f,interval,n=2):
    if n%2 != 0:
        raise ValueError
    mesh = np.linspace(interval[0],interval[1],n+1)
    count2 = int(n/2 - 1)
    coef = np.array([1]+[4,2]*count2+[4,1])
    h = (interval[1]-interval[0])/n
    def f_super(x):
        if abs(x) <= 1e-12:
            return f(1e-12)
        else:
            return f(x)
    fs = np.array([f_super(i) for i in mesh])
    result = h/3*(coef.dot(fs))
    return result
    
def doubleSimpson(f,lower,upper,interval,x_n,y_n):
    def f_sub(x):
        return compositeSimpson(partial(f,x),[lower(x),upper(x)],y_n)
    result = compositeSimpson(f_sub,interval,x_n)
    return result