# -*- coding: utf-8 -*-
"""
Created on Sat Sep 15 10:05:33 2018

@author: 49048
"""
from __future__ import division
from math import sqrt
import numpy as np
from matplotlib import pyplot as plt
#We've got to use recursion for ecaluating a_n
def recursion(n,x=2):
    #n is the recursion's depth
    if x < n:
        return sqrt(1+x*recursion(n,x+1))
    elif x == n:
        return sqrt(1+n)
    else:
        return None
        print 'Recursion times x has exceeded n!'
def a(n):
    #a(n) will return the value of the sequentce, \
    #by invoking the recursion and assigning n-depth to it
    if n==1:
        return 1
    elif n==2:
        return sqrt(1+2)
    elif n>=3:
        return recursion(n)
    else:
        return None
        print 'n must be a positive integer!'    
def main():
    for n in range(1,41):
        print '\hline n:%s\t & an:%s \\\\'%(n,a(n))    
    #guess a=3,plot ln(|an-a|)
    n = np.arange(1,41)
    logdeviation = np.array([np.log(abs(a(i)-3)) for i in n])
    line = np.array([3-(np.log(2))*i for i in n])
    figure = plt.figure()
    plt.plot(n,logdeviation,linewidth=2,label=r'$\ln(a_{n}-a)$')
    plt.plot(n,line,label=r'$y=3-(\ln2)n$')
    plt.legend()
    plt.grid()
    plt.show()
    #based on the figure, we can guess a_n=3+O(\beta_n),
    #where \beta_n=\frac{1}{2^n}
    
    
main()