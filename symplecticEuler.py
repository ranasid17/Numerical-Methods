#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import math 

# Symplectic Euler Method(s) for solving Hamiltonian Equations 
# Inputs
#   f : ODE (must be part of Hamiltonian)
#   g : ODE ( "     " )
#   q0: initial condition 1
#   p0: initial condition 2 
# Outputs
#   q : array containing num approx of Hamiltonian 
#   p : array "                                 " 

def symplecticEuler1(f, g, q0, p0, h):
    N = int(1/h) # calculate num steps from step size 
    q = np.zeros((N+1, np.size(q0))) # initialize output array 
    p = np.zeros((N+1, np.size(p0))) # initialize second output array 
    q[0] = q0
    p[0] = p0
    for n in range(N):
        q[n+1] = q[n] + h*f(p[n])
        p[n+1] = p[n] + h*g(q[n+1]) #modify second line to use q[n+1] term
    return q, p

# sE2 is identical to sE1 
# only vary in calculating q(n+1) or p(n+1) first 
def symplecticEuler2(f, g, q0, p0, h):
    N = int(1/h)
    q = np.zeros((N+1, np.size(q0)))
    p = np.zeros((N+1, np.size(p0)))
    q[0] = q0
    p[0] = p0
    for n in range(N):
        # Invert q, p lines 
        p[n+1] = p[n] + h*g(q[n]) 
        q[n+1] = q[n] + h*f(p[n+1]) #modify q to use p[n+1] term
    return q, p

# Sample Hamiltonian ODE (simple harmonic oscillator)
def f(p):
    return p
def g(q):
    return -q

# Sample Hamiltonian ODE (simple pendulum)
def func(p):
    return p
def func2(q):
    return -math.sin(q)