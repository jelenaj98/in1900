#Exercise E.41

import numpy as np
import matplotlib.pyplot as plt
from ODESolver import *

def SIR( u, t):             #function for solving differential eq. in the SIR model

#I = infected
#S = susceptibles
#R = recovered

    S,I,R = u
    nu = 0.1                #parameter in the ODE system
    dS = -beta*S*I          #S equation
    dI = beta*S*I - nu*I    #I equation
    dR = nu*I               #R equation
    return [dS,dI,dR]

def terminate(u,t,i):   #terminate condition- return true for termination if S + I + R is not sufficiently constant
    tol = 1e-6
    diff = (u[i,0]+u[i,1]+u[i,2])-(u[0,0]+u[0,1]+u[0,2])
    return abs(diff)>tol

def solve_SIR():
    solver = RungeKutta4(SIR)   #define object of RungeKutta4
    solver.set_initial_condition([1500,1,0])    #set initial conditions for RungeKutta4 object

    T = 60          #time of simulation [days]
    dt = 0.5        #minimum time interval(between two iterations)
    n = int(T/dt)   #number of iterations for T
    time = np.linspace(0,T,n+1)
    u, t = solver.solve(time, terminate)    #calculate u,t
    return u, t

def plot_SIR(t, u, filename): #plott the results
    S, I, R  = u[:,0], u[:,1], u[:,2]
    plt.plot(t,S,label = 'Susceptibles')
    plt.plot(t,I,label = 'infected')
    plt.plot(t,R,label = 'Recovered')
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('Spreading of a disease')
    plt.suptitle('Spreading of a disease by a SIR model')
    plt.savefig(filename)
    plt.show()

beta = 0.0005
u,t = solve_SIR()
plot_SIR(t,u,'SIR_0005.png')

beta = 0.0001
u,t = solve_SIR()
plot_SIR(t,u,'SIR_0001.png')

'''
run SIR.py

According to the equation for S:

S(t+dt)= S(t)-beta*S*I*dt

we can conclude that increasing beta will lead to increase number of infected and move from S to I category (S will rapidly decrease with increasing beta).
Contrary, with small beta, the number of infected will reduce (insignificant influence to S).
'''
