# CelestialDynamics
A Hamiltonian based planetary motion simulator.

![ten randomly generated plannets smashing into eachtoher](\media\collection.gif)

## Usage

1. Create instances of `celestial_object`.
```
earth = CelestialDynamics.celestial_object(earth_mass, [0,0,0], [0, 0, 0], colour='skyblue', radius=earth_radius)
moon = CelestialDynamics.celestial_object(moon_mass, [moon_distance, 0, 0], [0, moon_mass * moon_velocity, 0],
        colour='white', radius=moon_radius)
```
2. Define a `system` given by a list `celestial_object`'s
```
system = CelestialDynamics.stellar_system([earth, moon])
```
3. Integrate the system through time using `integrate`,
```
output = CelestialDynamics.integrate(system, t_eval, method='DOP853')
```
or, export the evolution directly to a video/gif using `animate_system_2D`.
```
CelestialDynamics.animate_system_2D(system, t_eval, 'earth_moon.gif')
```

## See Also
[1]  V. I. Arnold (1978). Mathematical Methods of Classical Mechanics. Springer New York, NY.

[2]  H. Goldstein, C. Poole, J. Safko (2002). Classical Mechanics (3rd Edition). Addison Wesley.

[3]  J. R. Taylor (2005). Classical Mechanics. University Science Books.

[4]  N. A. Lemos (2018). Analytical Mechanics. Cambridge University Press.

[5] SciPy documentation. https://docs.scipy.org/doc/scipy/
