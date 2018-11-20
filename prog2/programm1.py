# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 10:13:21 2018

@author: Alex
"""
from __future__ import division
from math import sqrt
import numpy as np
from matplotlib import pyplot as plt

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

def solve_for_absExtrema(polynomial,interval=[0,1]):
    candidates = [interval[0],interval[1]]
    #The possible absolute extrema, 
    #includiong the critical points and the end points of the interval
    critical_points = polynomial.get_derivative_solution()
    #Use the method of the polynomial to solve for its critical points
    if critical_points == 'No root!':
        pass#Do nothing
    elif type(critical_points)==float:
        if (critical_points>interval[0])&(critical_points<interval[1]):
            #check if the critical point is in the interval
           candidates.append(critical_points)
        else:
            pass#Do nothing
    else:
        for x in critical_points:
            if (x>interval[0])&(x<interval[1]):
                candidates.append(x)
            else:
                pass
    #Now we've decided the candidates, let's check which ones are extrema
    values = [polynomial.evaluate(x) for x in candidates]
    pairs = dict(zip(candidates,values))
    sort_result=sorted(pairs.items(),key=lambda x:x[1])
    return (sort_result[0],sort_result[-1])

def main():
    trials=(([-1,2],[-1,2,-1,1]),
            ([1,2],[1,-2,-1,1]),
            ([-2,1],[4,8,-4,-2]),
            ([-1,2],[1,0,1,-3]),
            ([-0.3,0.6],[1e-14,9,-3,0]),
            ([-1,2],[0,0,0,1.7]),
            ([0,3],[-3,9,-1e-14,0]),
            ([0,1],[0,-2,3,-1]))
    ntrials = 0
    for parameters in trials:
        print 'New trial\n Interval:%s \n Coefficients:%s'%(parameters[0],
                                                            parameters[1])
        polynomial = polynomial_deg3()
        polynomial.set_parameter(parameters[1])
        extrema = solve_for_absExtrema(polynomial,interval=parameters[0])
        print 'x_min:%.6f,p(x_min):%.6f'%(extrema[0][0],extrema[0][1])
        print 'x_max:%.6f,p(x_max):%.6f'%(extrema[1][0],extrema[1][1])
        print ''
        if ntrials in [2,6]:
            #If is the 3th or 7th trial, we will plot the figure.
            a = parameters[0][0]
            b = parameters[0][1]
            x = np.linspace(a,b,256,endpoint=True)
            y = np.array([polynomial.evaluate(i) for i in x])
            xmin = extrema[0][0]
            pmin = extrema[0][1]
            xmax = extrema[1][0]
            pmax = extrema[1][1]
            fig = plt.figure()
            plt.plot(x,y)
            plt.grid()
            #Now annotate the extrema
            plt.scatter(extrema[0][0],extrema[0][1],50,color='purple')
            plt.scatter(extrema[1][0],extrema[1][1],50,color='red')
            plt.annotate(r'$p_{min}=p(%.6f)=%.6f$'%(xmin,pmin),fontsize=10,
                         xy=(xmin,pmin),xytext=(10,30),
                         textcoords='offset points',
                         arrowprops=dict(arrowstyle='->',
                                         connectionstyle='arc3,rad=.3'))
            plt.annotate(r'$p_{max}=p(%.6f)=%.6f$'%(xmax,pmax),fontsize=10,
                         xy=(xmax,pmax),xytext=(10,-30),
                         textcoords='offset points',
                         arrowprops=dict(arrowstyle='->',
                                         connectionstyle='arc3,rad=-.1'))
            plt.show()
        ntrials = ntrials+1
        
        
main()

    
    
    