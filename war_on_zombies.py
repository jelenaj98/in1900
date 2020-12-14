import numpy as np
import matplotlib.pyplot as plt
from ODESolver import *
from math import exp
from SIZR import ProblemSIZR, SolverSIZR


#Test war on zombies
if __name__ == '__main__': # prevent calling thes methods from another modules e.g. using import
    T1, T2, T3 = 5, 10, 18 #strong atack phases
    T = [T1,T2,T3]
    beta = 0.03 #initil value of beta

    def calc_alpha(t):
        omega_small = 0.2*beta #small resistances against zombies from the humans throughout the simulation
        a = 50*omega_small
        sigma = 0.5
        omega = 0.0
        n = len(T)
        for i in range(n):  #calculating omega for according to strong attack phases
            omega += exp(-0.5*((t-T[i])/sigma)**2)

        omega *= a;
        omega += omega_small # add initial small resistance against zombies
        return omega

    problem = ProblemSIZR(
    beta = beta, S0 = 50, I0=0.0, Z0=3.0, R0=0, T=20,
    sigma = 0, deltaS=0, deltaI=0, ro=1,
    alpha=calc_alpha)   #alpha depending on T[T1,T2,T3] and calculates for each of them

    solver = SolverSIZR(problem, 0.5)   #use SIZR Model
    solver.solve()
    solver.plot('war on zobies', 'war on zobies')

'''
run war_on_zombies.py

The war on zombies modeled by the suggested omega(t) is sufficient
to save mankind :)
'''
