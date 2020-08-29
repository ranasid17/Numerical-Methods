#!/usr/bin/env python3

import numpy as np 

# Euler's Method Implementation for Partitioned ODEs
# Inputs
#   f : input function 1
#   g : input function 2 
#   q0, p0 : initial conditions for the respective functions 
#   N : number of steps 
# Outputs
#   q, p : numerical approximation to solution 
def euler(f, g, q0, p0, N):
    h = 1 / N # Calculate step size (inverse of num steps)
    q = np.zeros((N+1, np.size(q0))) # define array to hold ODE values
    p = np.zeros((N+1, np.size(p0))) # "  " 
    q[0] = q0 # store initial condition
    p[0] = p0 # "   " 
    # iteratively solve and store ODE
    for n in range(N):
        q[n+1] = q[n] + h*f(p[n])
        p[n+1] = p[n] + h*g(q[n])
    return q, p
