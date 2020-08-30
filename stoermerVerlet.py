#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np 
import math 

# Stoermer/Verlet Method for solving Hamiltonian Systems 
# Inputs
#   f : ODE (must be of Hamiltonian form)
#   g : ODE ( "                       " )
#   q0: initial condition 1
#   p0: initial condition 2 
#   h : step size
#   N : total number of steps 
# Outputs
#   q : array containing num approx of Hamiltonian 
#   p : array "     

def stoermerVerlet(f, g, q0, p0, h, N): 
    # initialize both output arrays 
    q = np.zeros((2*(N+1), np.size(q0)))
    p = np.zeros((2*(N+1), np.size(p0)))
    # store both initial conditions 
    q[0] = q0
    p[0] = p0
    # iteratively apply S/V method 
    for n in range(N):         
         pHalf = p[n] + ( 0.5 * h * g(q[n]) )
         q[n+1] = q[n] + ( h * f(pHalf) )
         p[n+1] = pHalf + ( 0.5 * h * g(q[n+1]) )
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