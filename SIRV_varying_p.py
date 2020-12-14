#Exercise E.44

import numpy as np
import matplotlib.pyplot as plt
from ODESolver import *
from SIRV import ProblemSIRV

class SolverSIRV:   #define class SolverSIRV
    def __init__(self, problem):
        self.problem = problem

    def solve(self, dt, method = RungeKutta4):  #solving the method
        self.solver = method(self.problem)
        self.solver.set_initial_condition(self.problem.getInitialCondition())
        u, self.t = self.solver.solve(self.problem.getTime(0,dt))
        self.S, self.I, self.R, self.V = u[:,0], u[:,1], u[:,2], u[:,3]

    def getMaxI(self):  #return the max value of I array
        return round(max(self.I))

    def plot(self): #plott the sulution
        plt.plot(self.t,self.S,label = 'Susceptibles')
        plt.plot(self.t,self.I,label = 'infected')
        plt.plot(self.t,self.R,label = 'Recovered')
        plt.plot(self.t,self.V,label = 'Vaccinated')
        plt.legend()
        plt.xlabel('Time')
        plt.ylabel('Spreading of a disease')
        plt.suptitle('Vaccination campaign in a SIR model')
        plt.savefig('SIRV_varying_p.png')
        plt.show()
        print('Maximum number of infected people is: %d' %(self.getMaxI()))


#Test SolverSIR
if __name__ == '__main__': # prevent calling thes methods from another modules e.g. using import
    problem = ProblemSIRV(beta = lambda t: 0.0005 if t<12 else 0.0001,  #beta is function of time
    nu = 0.1, S0 = 1500, I0=1, R0=0, T=60,V0=0,
    p = lambda t: 0.1 if t>=6 and t<=15 else 0)                       #p is function of time

    solver = SolverSIRV(problem)
    solver.solve(0.5)
    solver.plot()

'''
run SIRV_varying_p.py
Maximum number of infected people is: 388

Comment: Maximum number of infected people is 388 because vaccination lasts for only 10 days.
'''
