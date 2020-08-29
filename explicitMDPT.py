#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

# Explicit Midpoint (2nd order RK Method) for 1st order ODE
# Inputs
#   f :  input function
#   t0: initial time condition
#   y0: initial condition 
#   h : step size
# Outputs
#   y : array containing num approximations of each step (final index = soln)
def explicitMDPT(f, t0, y0, h): 
    # calc num steps from inverse of step size 
    N = int(1 / h) 
    # initialize time array 
    t = t0 + np.arange(N+1)*h 
    # initialize output array 
    y = np.zeros((N+1, np.size(y0))) 
    # store initial condition
    y[0] = y0 
    # Apply eMDPT iteratively
    for i in range(N): 
        # 1st intermediate stage and approximation
        xi1 = y[i]
        f1 = f(t[i], xi1)
        # 2nd intermediate stage and approximation
        xi2 = y[i] + h/2 * f1
        f2 = f(t[i+1] + h/2, xi2)
        # final approximation and storage
        y[i+1] = y[i] + h * f2
    # return array of num approximations 
    return y 

# Sample ODE for approximation
def model(t,y):
    dydt = y
    return dydt
# Sample IVP for approximation
t0 = 0
y0 = 1