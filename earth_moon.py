'''
Simulation of Earth, Moon System
'''
import numpy as np
import CelestialDynamics as cdp

earth_mass = cdp.kg_to_solar_mass(5.972e24)
earth_radius = 6371.0
moon_mass =  cdp.kg_to_solar_mass(7.34767309e22)
moon_distance = cdp.km_to_AU(384400)
moon_velocity = cdp.km_s_to_AU_yr(1.022)
moon_radius = 1737.0

earth = cdp.celestial_object(earth_mass, [0,0,0], [0, 0, 0], colour='skyblue', radius=earth_radius)
moon = cdp.celestial_object(moon_mass, [moon_distance, 0, 0],
                             [0, moon_mass * moon_velocity, 0], colour='white', radius=moon_radius)
system = cdp.stellar_system([earth, moon])

cdp.celestial_object
t_eval = np.linspace(0, 1/12, 1001)
# output = cdp.integrate(system, t_eval, method='DOP853')
cdp.animate_system_2D(system, t_eval, 'earth_moon.mp4')