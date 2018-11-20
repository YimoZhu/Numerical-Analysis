# -*- coding: utf-8 -*-
"""
Spyder Editor
Author: Yimo Zhu
This is for  Math 128A programming assignment #3
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

###############################################################################
def q1(a,b,x):
    #Take the input vectors a&b to use Recurrence System(3) to evaluate \phi at x
    #First lets check wether len(a)+1 = len(b). If not, the input is not valid
    if len(a)+1 != len(b):
        raise ValueError
    #Then if the input is valid:
    n = len(a)
    m = len(x)
    PHI = np.zeros((n+1,m))
    k = 0
    for x_k in x:
        #for each x we recursively compute \phi_0(x_k/) to \phi_n(x_k)
        phi = np.zeros(n+1)
        counter = 0
        while True:
            if counter == 0:
                phi[counter] = np.sqrt(1/2)
                counter = counter + 1
                continue
            elif counter == 1:
                phi[counter] = ((x_k-a[counter-1])*phi[counter-1])/np.sqrt(b[counter])
                counter = counter + 1
                continue
            elif counter <= n:
                phi[counter] = ((x_k-a[counter-1])*phi[counter-1]-np.sqrt(b[counter-1])*phi[counter-2])/np.sqrt(b[counter])
                counter = counter + 1
                continue
            else:
                break
        PHI[:,k] = phi
        k = k + 1
    return PHI


#---------------------------------------------------------------------
def q2(a,b):
    #Take the input vectors a&b to use Recurrence System(3) to get the coefficients of unit legendre polynomials up to degree n
    #First lets check wether len(a)+1 = len(b). If not, the input is not valid
    if len(a)+1 != len(b):
        raise ValueError
    #Then if the input is valid:    
    n = len(a)
    S = np.zeros((n+1,n+1))
    j = 0
    while True:
        if j == 0:
            S[j,j] = np.sqrt(1/b[j])
            j = j + 1
            continue
        elif j == 1:
            S[j,:] = (S[j-1,:]-np.append([0],a[j-1]*S[j-1,:-1]))/np.sqrt(b[j])
            j = j + 1
            continue
        elif j <= n:
            S[j,:] = (S[j-1,:]-np.append([0],a[j-1]*S[j-1,:-1])-np.append([0,0],np.sqrt(b[j-1])*S[j-2,:-2]))/np.sqrt(b[j])
            j = j + 1
            continue
        else:
            break
    return S


#---------------------------------------------------------------------
def q3(a,b):
    """Take the input vectors a&b to use olub-Welsch algorithm to calculate the abiscissas vector x & the weights vector w, each with length n"""
    #First lets check wether len(a)+1 = len(b). If not, the input is not valid
    if len(a)+1 != len(b):
        raise ValueError
    #Then if the input is valid:
    n = len(a)
    #Shape the coefficient matrix A
    A = np.zeros((n,n))
    row_counter = 0
    for a_i in a:
        A[row_counter,row_counter] = a_i
        if row_counter == 0:
            A[row_counter,row_counter+1] = np.sqrt(b[row_counter+1])
        elif row_counter == n-1:
            A[row_counter,row_counter-1] = np.sqrt(b[row_counter])
        else:
            A[row_counter,row_counter-1] = np.sqrt(b[row_counter])
            A[row_counter,row_counter+1] = np.sqrt(b[row_counter+1])
        #Finished filling the values in this row
        row_counter = row_counter + 1
        continue
    #Finished shaping the A. Now compute the eigenvalues
    x,Q = np.linalg.eig(A)
    w = b[0]*Q[0,:]**2
    return x,w


#---------------------------------------------------------------------
def q4(a,b):
    #Call the code from q2 to get the coefficients up to degree n
    S = q2(a,b)
    n = len(a)
    x = np.array(sorted(np.roots(S[n,:n+1])))
    Phi = np.zeros((n,n))
    for j in range(n):
        Phi[j,:] = np.polyval(S[j,:j+1],x)
    w = 1.0/np.diag(Phi.T.dot(Phi))
    return x,w
            

###############################################################################
#The main programms
    
"""Question # 1"""
#Now lets generate the 2 plots of \phi_n(x)
for n in [10,50]:
    print("Q1 Case n <-- %s"%n)
    a = np.zeros(n)
    b = np.append([2],np.arange(1,n+1)**2/(4*np.arange(1,n+1)**2-1))
    x = np.cos(np.linspace(-np.pi,0,10*n))
    phi_n = q1(a,b,x)[-1,:]
    fig = plt.figure(figsize=(10,8))
    plt.suptitle("Q1 Case n = %s"%n)
    plt.plot(x,phi_n)
    plt.grid(True)
    plt.show()    


#---------------------------------------------------------------------
"""Question # 2"""
n = 50
a = np.zeros(n)
b = np.append([2],np.arange(1,n+1)**2/(4*np.arange(1,n+1)**2-1))
S = q2(a,b)
x = np.cos(np.linspace(-np.pi,0,500))
phi10 = S[10,:11]
phi50 = S[50,:51]
#First the plot where n = 10
fig = plt.figure(figsize=(10,8))
plt.suptitle("Q2 Case n = 10")
y10 = np.polyval(phi10,x)
plt.plot(x,y10)
plt.grid(True)
fig = plt.figure(figsize=(10,8))
plt.suptitle("Q2 Case n = 50")
y50 = np.polyval(phi50,x)
plt.plot(x,y50)
plt.grid(True)
plt.show()


#---------------------------------------------------------------------
"""Question # 3"""
n = 10
E = np.zeros(2*n+1)
a = np.zeros(n)
b = np.append([2],np.arange(1,n+1)**2/(4*np.arange(1,n+1)**2-1))
x,w = q3(a,b)
for k in np.arange(0,2*n+1):
    E[k] = abs(w.dot(np.cos(k*np.arccos(x))) - (1+(-1)**k)/(1-k**2 +1.0e-18))
xw = pd.DataFrame({'abscissas':x,'weight':w})
print(xw)
xw.to_csv('Q3xwCase10.csv')
print(E)
pd.DataFrame({'E':E}).to_csv('Q3ECase10.csv')

n = 40
E = np.zeros(2*n+1)
a = np.zeros(n)
b = np.append([2],np.arange(1,n+1)**2/(4*np.arange(1,n+1)**2-1))
x,w = q3(a,b)
for k in np.arange(0,2*n+1):
    E[k] = abs(w.dot(np.cos(k*np.arccos(x))) - (1+(-1)**k)/(1-k**2 +1.0e-18))
xw = pd.DataFrame({'abscissas':x,'weight':w})
print(xw)
xw.to_csv('Q3xwCase40.csv')
print(E)
pd.DataFrame({'E':E}).to_csv('Q3ECase40.csv')
print('Norm:%s'%np.linalg.norm(E[:80]))


#---------------------------------------------------------------------
"""Question # 4"""
n = 10
E = np.zeros(2*n+1)
a = np.zeros(n)
b = np.append([2],np.arange(1,n+1)**2/(4*np.arange(1,n+1)**2-1))
x,w = q4(a,b)
for k in np.arange(0,2*n+1):
    E[k] = abs(w.dot(np.cos(k*np.arccos(x))) - (1+(-1)**k)/(1-k**2 +1.0e-18))
xw = pd.DataFrame({'abscissas':x,'weight':w})
print(xw)
xw.to_csv('Q4xwCase10.csv')
print(E)
pd.DataFrame({'E':E}).to_csv('Q4ECase10.csv')

n = 40
E = np.zeros(2*n+1)
a = np.zeros(n)
b = np.append([2],np.arange(1,n+1)**2/(4*np.arange(1,n+1)**2-1))
x,w = q4(a,b)
for k in np.arange(0,2*n+1):
    E[k] = abs(w.dot(np.cos(k*np.arccos(x))) - (1+(-1)**k)/(1-k**2 +1.0e-18))
xw = pd.DataFrame({'abscissas':x,'weight':w})
print(xw)
xw.to_csv('Q4xwCase40.csv')
print(E)
pd.DataFrame({'E':E}).to_csv('Q4ECase40.csv')
print('Norm:%s'%np.linalg.norm(E[:80]))