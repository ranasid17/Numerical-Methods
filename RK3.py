#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np 

# Cannonical 3rd order RK method for solving ODEs
# Inputs
#   f : ODE function
#   t0: initial time condition
#   y0: initial condition
#   h : step size 
# Outputs
#   y : array of num approximations for each step 
def rk3(f, t0, y0, h): 
    # define num steps from inverse of step size
    N = int(1/h)
    # initialize time array 
    t = t0 + np.arange(N+1)*h 
    # initialize output array
    y = np.zeros((N+1, np.size(y0))) 
    # store IC 
    y[0] = y0 
    # apply RK3 iteratively  
    for i in range(N): 
        # 1st intermediate stage + evaluation
        xi1 = y[i]
        f1 = f(t[i], xi1)
        # 2nd intermediate stage + evaluation
        xi2 = y[i] + (h/2 * f1)
        f2 = f(t[i]+ h/2, xi2)
        # 3rd intermediate stage + evaluation
        xi3 = y[i] + h * (-f1 + 2*f2)
        f3 = f(t[i+1] + h, xi3)
        # final approximation + storage
        y[i+1] = y[i] + (h/6) * (f1 + 4*f2 + f3) 
    return y 
        
# Sample ODE for approximation
def model(t,y):
    dydt = y
    return dydt
# Sample IVP for approximation
t0 = 0
y0 = 1
