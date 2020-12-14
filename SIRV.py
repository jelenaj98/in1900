#Exercise E.43

import numpy as np
import matplotlib.pyplot as plt
from ODESolver import *
from SIR_class import ProblemSIR

class ProblemSIRV(ProblemSIR):  #define class for the SIRV model as a subclass of ProblemSIR
    def __init__(self, nu, beta,  S0, I0, R0, T, V0, p):
        ProblemSIR.__init__(self,nu, beta, S0, I0, R0, T)   #calling __init__ from Superclass
        self.V0= V0     #initialization of vaccination (belongs to ProblemSIRV)

        # wrap fraction as a function of time or constant
        if isinstance(p,(float,int)):
            self.p= lambda t: p
        elif callable(p):
            self.p= p


    def __call__(self, u, t):   #calculating equation for SIRV model
        S,I,R,V = u
        return [-self.beta(t)*S*I - S*self.p(t),
        self.beta(t)*S*I - self.nu(t)*I,
        self.nu(t)*I,
        self.p(t)*S]

    def terminate(self,u,t,i):  #terminate condition- return true for termination if S + I + R + V is not sufficiently constant
        tol = 1e-6
        diff = (u[i,0]+u[i,1]+u[i,2]+u[i,3])-(u[0,0]+u[0,1]+u[0,2]+u[0,3])
        return abs(diff)>tol

    def solve(self, dt):
        self.solver = RungeKutta4(self) #Solve the RungeKutta4 equation
        self.solver.set_initial_condition([1500,1,0,0]) #set initial values for S,I,R and V
        u, self.t = self.solver.solve(self.getTime(0,dt))
        self.S, self.I, self.R, self.V = u[:,0], u[:,1], u[:,2], u[:,3] #create arrays for S,I, R and V

    def getInitialCondition(self):
        return [self.S0, self.I0, self.R0, self.V0]

    def plot(self): #plott the results
        plt.plot(self.t,self.S,label = 'Susceptibles')
        plt.plot(self.t,self.I,label = 'infected')
        plt.plot(self.t,self.R,label = 'Recovered')
        plt.plot(self.t,self.V,label = 'Vaccinated')
        plt.legend()
        plt.xlabel('Time')
        plt.ylabel('Spreading of a disease')
        plt.suptitle('Vaccination in a SIRV model')
        plt.savefig('SIRV.png')
        plt.show()
        print('Maximum number of infected people is: %d' %(round(max(self.I))))


if __name__ == '__main__': # prevent calling thes methods from another modules e.g. using import
    problem = ProblemSIRV(beta = lambda t: 0.0005 if t<12 else 0.0001,
    nu = 0.1, S0 = 1500, I0=1, R0=0, T=60, V0 = 0, p=0.1)
    problem.solve(dt=0.5)
    problem.plot()

'''
run SIRV.py
Maximum number of infected people is: 50

Comment: Vaccination effect will lead to decreasing number of infected. Max number of infected is 50
which comparesing with previous model (732) is much lower!!
'''
