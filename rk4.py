#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

# "The" cannonical RK method (RK4)
# Inputs
#   f : ODE
#   t0: initial time condition
#   y0: initial condition
#   h : step size
# Outputs
#   y : array holding approximations
def rk4(f, t0, y0, h): 
    # define num steps from inverse step size 
    N = int(1/h)
    # initialize time array 
    t = t0 + np.arange(N+1)*h 
    # initialize output array
    y = np.zeros((N+1, np.size(y0)))
    # store initial condition
    y[0] = y0    
    # iterate 
    for i in range(N): 
        # 1st intermediate stage and evaluation
        xi1 = y[i]
        f1 = f(t[i], xi1)
        # 2nd "    "
        xi2 = y[i] + (h * 0.5 * f1) 
        f2 = f(t[i] + h * 0.5, xi2)
        # 3rd "    "
        xi3 = y[i] + (h * 0.5 * f2)
        f3 = f(t[i] + h * 0.5, xi3)
        # 4th "    "
        xi4 = y[i] + h * f3
        f4 = f(t[i] + h, xi4)
        # final approximation and storage
        y[i+1] = y[i] + (h * 1/6) *(f1 + 2*f2 + 2*f3 + f4)                
   
    return y

# Sample ODE for approximation
def model(t,y):
    dydt = y
    return dydt
# Sample IVP for approximation
t0 = 0
y0 = 1
