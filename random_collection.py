'''
Simulation of 10 Planets with Randomised Initial Conditions
'''
import numpy as np
import CelestialDynamics as cdp

n = 10
stars = []

for i in range(n):
    mass_i = np.random.uniform(1, 100)
    pos_i = np.random.uniform(-5, 5, 3)
    momenta_i = np.random.uniform(-20, 20, 3) * mass_i
    stars.append(cdp.celestial_object(mass_i, pos_i, momenta_i, mass_i)) 

system = cdp.stellar_system(stars)
t_eval = np.linspace(0, 1.2, 701)

cdp.animate_system_2D(system, t_eval, 'collection.gif')


