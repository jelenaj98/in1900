#Exercise E.42

import numpy as np
import matplotlib.pyplot as plt
from ODESolver import *

class ProblemSIR:   #define class for the SIR model
    def __init__(self, nu, beta, S0, I0, R0, T):
        # wrap nu as a function of time or constant
        if isinstance(nu,(float,int)):
            self.nu= lambda t: nu
        elif callable(nu):
            self.nu= nu
        # wrap beta as a function of time or constant
        if isinstance(beta,(float,int)):
            self.beta= lambda t: beta
        elif callable(beta):
            self.beta= beta

        self.S0, self.I0, self.R0, self.T= S0, I0, R0, T    #set initial values

    def getTime(self,t0,dt):     #calculate and return evently spaced numbers or a specificed interval T
        n = int(self.T/dt)
        return np.linspace(t0,self.T,n+1)

    def __call__(self, u, t):   #calculating equation for SIR model
        S,I,R = u
        return [-self.beta(t)*S*I,
        self.beta(t)*S*I - self.nu(t)*I,
        self.nu(t)*I]

    def getInitialCondition(self):
        return [self.S0, self.I0, self.R0]

    def terminate(self,u,t,i):  #terminate condition- return true for termination if S + I + R is not sufficiently constant
        tol = 1e-6
        diff = (u[i,0]+u[i,1]+u[i,2])-(u[0,0]+u[0,1]+u[0,2])
        return abs(diff)>tol

class SolverSIR:    #definition of class SolverSIR
    def __init__(self, problem, dt): #problem: object of class ProblemSIR
        self.problem, self.dt = problem, dt

    def solve(self, method = RungeKutta4):  #Solve the RungeKutta4 equation
        self.solver = method(self.problem)
        self.solver.set_initial_condition(self.problem.getInitialCondition())
        u, self.t = self.solver.solve(self.problem.getTime(0,self.dt))
        self.S, self.I, self.R = u[:,0], u[:,1], u[:,2] #create arrays for S,I and R

    def plot(self):
        plt.plot(self.t,self.S,label = 'Susceptibles')
        plt.plot(self.t,self.I,label = 'infected')
        plt.plot(self.t,self.R,label = 'Recovered')
        plt.legend()
        plt.xlabel('Time')
        plt.ylabel('Spreading of a disease')
        plt.suptitle('Spreading of a disease by a SIR model')
        plt.savefig('SIR_class.png')
        plt.show()
        print('Maximum number of infected people is: %d' %(round(max(self.I))))

# initialization
problem = ProblemSIR(beta = lambda t: 0.0005 if t<12 else 0.0001,   # beta is function of t(beta will reduce after 12 days)
nu = 0.1, S0 = 1500, I0=1, R0=0, T=60)


#Test SolverSIR
if __name__ == '__main__': # prevent calling thes methods from another modules (stuff only to run when not called via 'import' here)
    solverSIR = SolverSIR(problem, 0.5)
    solverSIR.solve()
    solverSIR.plot()

'''
run SIR_class.py

Maximum number of infected people is: 732

Comment: During beta = 0.0005 the number of infected is growing and reaches the maximum,
but people will soon be aware of infenction and will take precautions (beta = 0.0001) which
leads to decreasing number of infected.
'''
