__version__ = "1.0"
__author__ = "Yimo Zhu"
docker = "A package for numerical analysis algorithms. As some algorithms put in some big packages are hard to find, \
          this one is specific and for some daily use."

import numpy as np
from numpy import linalg
import inspect
from time import ctime
from matplotlib import pyplot as plt

#Configuration
settings = type("pyNA_settingsObject",(),{"print_on":True,"figure_on":False})

def __output(msg):
    if settings.print_on == True:
        currentFrame = inspect.currentframe()
        calFrame = inspect.getouterframes(currentFrame,2)
        print("[%s] @%s : %s"%(ctime(),calFrame[1][3],msg))
    else:
        pass

def __flag_plot(func):
    if settings.figure_on == True:
        func()
    

def set_fields(fieldName,value):
    setattr(settings,fieldName,value)

def rungeKutta4(f,initial,a,b,N,mode="scaler"):
    #Parameter Notation:
    #f: the function object.
    #initial: Initial values of the IVP
    #a,b: the left and right boundary of the interval
    #N: number of nodes.
    #mode: Scaler or Vector
    if mode == "scaler":
        pass
    elif mode == "vector":
        d = len(initial)
        h = (b-a)/N
        t = np.linspace(a,b,N+1)
        w = np.empty((N+1,d))
        w[0,:] = initial
        for j in np.arange(N):
            k1 = f(t[j] , w[j,:])
            k2 = f(t[j]+0.5*h , w[j,:]+0.5*h*k1)
            k3 = f(t[j]+0.5*h , w[j,:]+0.5*h*k2)
            k4 = f(t[j]+h , w[j,:]+h*k3)
            w[j+1,:] = w[j,:] + h/6 * (k1 + 2*k2 + 2*k3 + k4)
            __output("Step # %s, result %s"%(j+1,w[j+1,:]))
        return w

def euler(f,initial,a,b,N,mode="scaler"):
    #Euler's method
    #The specification of input parameters are much the same as rungeKutta4
    if mode == "scaler":
        h = (b-a)/N
        t = np.linspace(a,b,N+1)
        w = np.empty(N+1)
        w[0] = initial
        for j in np.range(N):
            w[j+1] = w[j] + h*f(t[j],w[j])
            __output("Step # %s, result %s"%(j+1,w[j+1]))
        return w
    elif mode == "vector":
        pass
