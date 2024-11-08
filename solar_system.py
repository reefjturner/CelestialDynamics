'''
Simulation of Sun, Earth, Jupiter, Saturn System (introduce a rogue black hole too if you like)
'''
import numpy as np
import CelestialDynamics as cdp

earth_mass = cdp.kg_to_solar_mass(5.972e24)
earth_radius = 1
earth_period = 1
earth_velocity = cdp.km_s_to_AU_yr(29)

jupiter_mass = cdp.kg_to_solar_mass(1.898e27)
jupiter_radius = 5.2
jupiter_period = 11.86
jupiter_velocity = cdp.km_s_to_AU_yr(12.6)

saturn_mass = cdp.kg_to_solar_mass(5.6834e26)
saturn_radius = 9.5826
saturn_period = 29.4475
saturn_velocity = cdp.km_s_to_AU_yr(9.87)

sun = cdp.celestial_object(1, [0,0,0], [0,0,0], colour='yellow', radius=696340)
earth = cdp.celestial_object(earth_mass, [earth_radius,0,0], [0, earth_velocity * earth_mass, 0], colour='skyblue', radius=6371)
jupiter = cdp.celestial_object(jupiter_mass, [jupiter_radius, 0, 0], [0, jupiter_velocity * jupiter_mass, 0], colour='orange', radius=69911)
saturn = cdp.celestial_object(saturn_mass, [saturn_radius, 0, 0], [0, saturn_velocity * saturn_mass, 0], colour='green', radius=58232)
black_hole = cdp.celestial_object(10, [-15, -15, 0], [70, 130, 0], radius = 696340*3)

# system = cdp.stellar_system([sun, earth, jupiter, saturn])
system = cdp.stellar_system([sun, earth, jupiter, saturn, black_hole])

t_eval = np.linspace(0, 3, 301)
output = cdp.integrate(system, t_eval, method='DOP853')

cdp.animate_system_2D(system, t_eval, 'sun_earth_jupiter_ASTROID_LOOK_OUT.gif')