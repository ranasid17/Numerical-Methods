#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

# Heun's Method (2nd order RKM/Exp Trap) for solving ODEs
# Inputs
#   f : ODE
#   t0: initial time condition
#   h : step size 
#   y0: array 
# Output
#   y : array of step-wise approximations

def HeunsMethod(f, t0, y0, h):
    # calculate num steps from inverse step size 
    N = int(1/h)
    # initialize time array 
    t = t0 + np.arange(N+1)*h
    # initialize output array 
    y = np.zeros((N+1, np.size(y0)))
    # store initial value 
    y[0] = y0
    # iteratively apply method 
    for n in range(N):
        # first intermediate step  + evaluation 
        x1 = y[n]
        f1 = f(t[n], x1)
        # second intermediate step + evaluation 
        x2 = y[n] + h*f1
        f2 = f(t[n+1], x2)
        # store final approx 
        y[n+1] = y[n] + 0.5*h*(f1 + f2)
    return y

# Sample ODE for approximation
def model(t,y):
    dydt = y
    return dydt

# Sample IVP for approximation
t0 = 0
y0 = 1