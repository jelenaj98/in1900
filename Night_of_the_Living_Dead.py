#Exercise E.47

import numpy as np
import matplotlib.pyplot as plt
from ODESolver import *
from SIZR import ProblemSIZR, SolverSIZR


'''
sigma: the number of new humans brought into the zombified area per unit time
beta: the probability that a theoretically possible human-zombie pair actually meets physically, during a unit time interval, with the result that the human is infected
deltaS: the probability that a susceptible human is killed or dies, in a unit time interval
deltaI : the probability that an infected human is killed or dies, in a unit time interval
ro: the probability that an infected human is turned into a zombie, during a unit time interval
alpha: the probability that, during a unit time interval, a theoretically possible human-zombie pair fights and the human kills the zombie
'''

class piecewise:    #define piecewise class(for chosing value depending on time)
    def __init__(self,time,values):
        if len(time) != len(values):
            print("Wrong arguments to piecewise function")
        self.time = time
        self.values = values

    def __call__(self,t):
        for t_, v_ in zip(self.time,self.values):
            if t < t_:
                return v_
        return 0

#Test Night of the Living dead
if __name__ == '__main__': # prevent calling thes methods from another modules e.g. using import

    T1, T2, T3 = 4, 28, 33 #zombification phases

    problem = ProblemSIZR(
    beta = piecewise([T1,T2,T3],[0.03, 0.0012, 0.0]),
    S0 = 60, I0=0, Z0=1, R0=0, T=33,
    sigma = piecewise([T1,T2,T3],[20.0, 2.0, 0.0]),
    deltaS=piecewise([T1,T2,T3],[0.0, 0.0, 0.0067]),
    deltaI=piecewise([T1,T2,T3],[0.0, 0.014, 0.0]),
    ro=1,
    alpha=piecewise([T1,T2,T3],[0.0, 0.0016, 0.006]))

    solver = SolverSIZR(problem, 0.5)
    solver.solve()
    solver.plot('Night of the Living Dead', 'Night of the Living Dead')

'''
run Night_of_the_Living_Dead.py
'''
