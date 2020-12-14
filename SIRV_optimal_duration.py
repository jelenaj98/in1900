#Exercise E.45

import numpy as np
import matplotlib.pyplot as plt
from ODESolver import *
from SIRV_varying_p import SolverSIRV
from SIRV import ProblemSIRV


n = 31 #Vt = [0,31]
Vt,It = [],[] #initialise empty lists
negligibleEffect = 0.1 #negligible difference betven infected people (current - previous day)

#Compute the maximum number of infected people, maxt I(t), as a function of VT [0, 31]
def computeOprimalVaccinationPeriod():
    for i in range(n+1):
        Vt.append(i) # append to Vt before usage in a following row

        problem = ProblemSIRV(beta = lambda t: 0.0005 if t<12 else 0.0001,  #beta is function of time
        p = lambda t: 0.1 if t>=6 and t<=6+Vt[i] else 0.0,  #p is function of time and Vt
        nu = 0.1, S0 = 1500, I0=1, R0=0,V0=0, T=60)

        solver = SolverSIRV(problem)
        solver.solve(dt=0.5)

        It.append(solver.getMaxI())


def plotSolution():
    plt.plot(Vt,It,label = 'infected')
    plt.legend()
    plt.xlabel('Vaccination Period')
    plt.ylabel('Spreading of a disease')
    plt.suptitle('Optimal vaccination period')
    plt.savefig('SIRV_optimal_duration.png')
    plt.show()


# calculate the smallest vaccination period VT such that increasing VT
#has negligible effect on the maximum number of infected people.
def calculateSmallestVaccinationPeriod():
    for i in range(1,n+1):
        if abs(It[i]-It[i-1]) <= negligibleEffect:
            return i-1;

#Test
computeOprimalVaccinationPeriod()
calculateSmallestVaccinationPeriod()
plotSolution()
print('The Smallest vaccination period: %d days' %calculateSmallestVaccinationPeriod())

'''
run SIRV_optimal_duration.py

The Smallest vaccination period: 6 days

'''
