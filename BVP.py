#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.linalg import solveh_banded
import matplotlib.pyplot as plt

# Two point BVP integraters
# solveBVP1 solves PDE of form, -u"(x) = f(x), 0 < x < 1, u(0) = u(1)= 0
# solveBVP2 solves PDE of form, -u"(x) + u(x) = f(x), 0 < x < 1, u(0) = u(1) = 0
# Sample PDEs provided at end 

# Inputs 
#   f : ODE, must be of specified form above
#   N : numer of steps 
def solveBVP1(f,N): 
    # calc step size 
    h = 1/N
    # set boundary conditions
    y0 = 0
    yN = 0
    
    # define nodes
    x = np.linspace(0, 1, N)
    # define internal nodes
    x_int = x[1:N-1]
    
    # create LHS (banded) matrix 
    n_mat = N - 2
    Ah = np.zeros((2, n_mat))
    # Fill banded matrix Ab 
    Ah[0][:] = -1 
    Ah[1][:] = 2 
    # multiply Ah by 1/h^2
    Ah = Ah * (1/h**2)
        
    # create RHS matrix 
    b = np.zeros((N-2,1))
    # fill RHS with f(x_int)
    for i in range(n_mat):
        b[i] = f(x_int[i])
 
    # solve system 
    y_int = solveh_banded(Ah,b)
    output = np.append(y0, y_int)
    output = np.append(output, yN)
    
    # plot solns 
    plt.plot(x,output)
    plt.title("Fin Diff Approximation for " + f.__name__)
    plt.title("N =" + str(N), loc = 'left')

    return output 

# Inputs 
#   f : ODE, must be of specified form above
#   N : numer of steps 
def solveBVP2(f, N):
    # calc step size 
    h = 1/N
    # set boundary conditions
    y0 = 0
    yN = 0
    # define nodes
    x = np.linspace(0, 1, N)
    # define internal nodes
    x_int = x[1:N-1]
    
    # create LHS (banded) matrix 
    n_mat = N - 2
    Ah = np.zeros((2, n_mat))
    # Fill banded matrix Ab 
    Ah[0][:] = -1 
    Ah[1][:] = 2 - (h**2) 
    
    # create RHS matrix 
    b = np.zeros((N-2,1))
    # fill RHS with f(x_int)
    for i in range(n_mat):
        b[i] = f(x_int[i]) * h**2
        
    # solve system 
    y_int = solveh_banded(Ah,b)
    output = np.append(y0, y_int)
    output = np.append(output, yN)
    
    # plot solns 
    plt.plot(x,output)
    plt.title("Fin Diff Approximation for " + f.__name__)
    plt.title("N =" + str(N), loc = 'left')

    return output 

# Sample ODEs
def f1(x):
    return 1

def f2(x): 
    return x**2
