#Exercise E.46

import numpy as np
import matplotlib.pyplot as plt
from ODESolver import *
'''
sigma: the number of new humans brought into the zombified area per unit time
beta: the probability that a theoretically possible human-zombie pair actually meets physically, during a unit time interval, with the result that the human is infected
deltaS: the probability that a susceptible human is killed or dies, in a unit time interval
deltaI : the probability that an infected human is killed or dies, in a unit time interval
ro: the probability that an infected human is turned into a zombie, during a unit time interval
alpha: the probability that, during a unit time interval, a theoretically possible human-zombie pair fights and the human kills the zombie
'''

class ProblemSIZR:
    def __init__(self, beta, S0, I0, Z0, R0, T, sigma, deltaS,deltaI,ro,alpha):

        # wrap beta as a function of time or constant
        if isinstance(beta,(float,int)):
            self.beta= lambda t: beta
        elif callable(beta):
            self.beta= beta

        # wrap sigma as a function of time or constant
        if isinstance(sigma,(float,int)):
            self.sigma= lambda t: sigma
        elif callable(sigma):
            self.sigma = sigma

        # wrap deltaS as a function of time or constant
        if isinstance(deltaS,(float,int)):
            self.deltaS= lambda t: deltaS
        elif callable(deltaS):
            self.deltaS = deltaS

        # wrap deltaI as a function of time or constant
        if isinstance(deltaI,(float,int)):
            self.deltaI= lambda t: deltaI
        elif callable(deltaI):
            self.deltaI = deltaI

        # wrap ro as a function of time or constant
        if isinstance(ro,(float,int)):
            self.ro= lambda t: ro
        elif callable(ro):
            self.ro = ro

        # wrap alpha as a function of time or constant
        if isinstance(alpha,(float,int)):
            self.alpha= lambda t: alpha
        elif callable(alpha):
            self.alpha = alpha

        self.S0, self.I0, self.Z0, self.R0, self.T= S0, I0, Z0, R0, T #initialisation

    def getTime(self,t0,dt):        #calculate and return evently spaced numbers or a specificed interval T
        n = int(self.T/dt)
        return np.linspace(t0,self.T,n+1)

    def __call__(self, u, t):   #calculating equation for SIZR model
        S,I,Z,R = u
        return [self.sigma(t) - self.beta(t)*S*Z - self.deltaS(t)*S,    #S'
        self.beta(t)*S*Z - self.ro(t)*I - self.deltaI(t)*I,         #I'
        self.ro(t)*I - self.alpha(t)*S*Z,                           #Z'
        self.deltaS(t)*S + self.deltaI(t)*I + self.alpha(t)*S*Z]   #R'

    def getInitialCondition(self):
        return [self.S0, self.I0, self.Z0, self.R0]

    def terminate(self,u,t,i):  #terminate condition- return true for termination if S + I + Z + R is not sufficiently constant
        tol = 1e-6
        diff = (u[i,0]+u[i,1]+u[i,2]+u[i,3])-(u[0,0]+u[0,1]+u[0,2]+u[i,3])
        return abs(diff)>tol

class SolverSIZR:   #define class for solving SIZR model
    def __init__(self, problem, dt):
        self.problem, self.dt = problem, dt

    def solve(self, method = RungeKutta4):
        self.solver = method(self.problem)
        self.solver.set_initial_condition(self.problem.getInitialCondition())
        u, self.t = self.solver.solve(self.problem.getTime(0,self.dt))
        self.S, self.I, self.Z, self.R = u[:,0], u[:,1], u[:,2], u[:,3]

    def plot(self,title, filename):   #plotting the results
        plt.plot(self.t,self.S,label = 'Susceptibles humans')
        plt.plot(self.t,self.I,label = 'infected humans')
        plt.plot(self.t,self.Z,label = 'Zombies')
        plt.plot(self.t,self.R,label = 'Removed individuals')
        plt.legend()
        plt.xlabel('Time (h)')
        plt.ylabel('Zombification')
        plt.suptitle(title)
        plt.savefig(filename+'.png')
        plt.show()

#Test SolverSIZR
if __name__ == '__main__': # prevent calling thes methods from another modules e.g. using import
    problem = ProblemSIZR(beta = 0.0012,
    S0 = 10, I0=0, Z0=100, R0=0, T=24,
    sigma=2, deltaS=0.0, deltaI=0.014, ro=1, alpha=0.0016)

    solver = SolverSIZR(problem, 0.5)
    solver.solve()
    solver.plot('Simulate human-zombie interaction', 'SIZR')

'''
run SIZR.py
'''
