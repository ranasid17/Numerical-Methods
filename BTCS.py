#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pylab import *
from mpl_toolkits.mplot3d import Axes3D

# Backwards Time, Center Space method for solving system of PDEs
# BTCS applies Rev Euler's Method in time, centered finite diff in space

# Inputs
#   M : number of space steps
#   N : number of time steps 
# Output
#   Plot of numerical solution to PDE 

def plotBTCS(M,N): 
    # define space and time mesh 
    x = np.linspace(0,1,M+1)    
    t = np.linspace(0,1,N+1)    
    X, T = np.meshgrid(x,t,indexing='ij')    
    # define space and time steps   
    h = 1/M    
    k = 1/N    
    # define initial conditions for t0
    u = np.zeros((M+1,N+1))    
    u[:,0] = sin(pi*x)**2    
    # initialize finite difference matrix in space    
    Ah = (2*np.eye(M-1) - np.eye(M-1,k=1) - np.eye(M-1,k=-1))/h**2    
    # iteratively apply method 
    for n in range(N): 
       u[1:M,n+1] = u[1:M,n] - k*Ah @ u[1:M,n+1]    
    # plot solution    
    fig = plt.figure()    
    ax = Axes3D(fig)    
    ax.plot_surface(X,T,u,rstride=int(ceil(M/10)),cstride=int(ceil(N/10)))
    plt.title("M =" + str(M), loc = 'left')
    plt.title("N =" + str(N), loc = 'right')