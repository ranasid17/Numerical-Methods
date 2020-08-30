#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pylab import *
from mpl_toolkits.mplot3d import Axes3D

# Forwards Time, Center Space method for solving system of PDEs
# FTCS applies Euler's Method in time, centered second order fin diff in space
# Inputs
#   M : number of space steps
#   N : number of time steps 
# Output
#   Plot of numerical solution to PDE 

def FTCS(M,N):    
    # set up x and t grids    
    x = linspace(0,1,M+1)    
    t = linspace(0,1,N+1)    
    X, T = meshgrid(x,t,indexing='ij')    
    # set space step size h and time step size k    
    h = 1/M    
    k = 1/N    
    # set initial conditions at t=0    
    u = zeros((M+1,N+1))    
    u[:,0] = sin(pi*x)**2    
    # finite difference matrix in space    
    Ah = (2*eye(M-1) - eye(M-1,k=1) - eye(M-1,k=-1))/h**2    
    # iterate forward from t=0 to t=1    
    for n in range(N):        
        u[1:M,n+1] = u[1:M,n] - k*Ah @ u[1:M,n]    
    # plot solution    
    fig = figure()    
    ax = Axes3D(fig)    
    ax.plot_surface(X,T,u,rstride=int(ceil(M/10)),cstride=int(ceil(N/10)))
    plt.title("M =" + str(M), loc = 'left')
    plt.title("N =" + str(N), loc = 'right')

    
