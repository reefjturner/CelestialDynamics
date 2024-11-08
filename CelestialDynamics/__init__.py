'''
Documentation!
'''
import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as plt
import matplotlib.animation as animation

G = 39.478 # gravitational constant in AU, yr, Solar Mass

def km_s_to_AU_yr(v):
    '''
    Converts from km/sec to AU/year
    '''
    return v * 0.210805

def km_to_AU(s):
    return s * 6.68459e-9

def AU_to_km(s):
    return s / 6.68459e-9

def kg_to_solar_mass(m):
    return m * 5.02785e-31

class celestial_object:
    '''
    Stuff!

    Attributes:
        mass (float): Mass of the object.

        position (array_like, default = np.array([0.0,0.0,0.0])): 1-D arrays of length 3 which represents the position of the object in AU.

        dot_position (array_like, default = np.array([0.0,0.0,0.0])): 1-D arrays of length 3 which represents the time derivative of the position of the object in AU / year.

        momentum (array_like, default = np.array([0.0,0.0,0.0])): 1-D arrays of length 3 which represents the momentum of the object in solar masses AU / year.

        dot_momentum (array_like, default = np.array([0.0,0.0,0.0])): 1-D arrays of length 3 which represents the time derivative of the momentum of the object in solar masses AU / year^2.

        radius (float, default = 1.0): Radius of object in km.

        colour (str, default = 'white'): Colour of object.

        temp (float, default = 300.0): Tempreature of object in K.

    Methods:
        To be done still :P
    '''
    def __init__(self, mass: float, position=np.array([0.0,0.0,0.0]), momentum=np.array([0.0,0.0,0.0]), radius=1.0, colour='white', temp=300.0):
        '''
        Bruh!

        Parameters:
            mass (float): Mass of the object.

            position (array_like, default = np.array([0.0,0.0,0.0])): 1-D arrays of length 3 which represents the initial position of the object in AU.

            momentum (array_like, default = np.array([0.0,0.0,0.0])): 1-D arrays of length 3 which represents the initial momentum of the object in solar masses AU / year.

            radius (float, default = 1.0): Radius of object in km.

            colour (str, default = 'white'): Colour of object.

            temp (float, default = 300.0): Tempreature of object in K.

        Examples:
        >>> earth_mass = CelestialDynamics.kg_to_solar_mass(5.972e24)
        >>> earth_radius = 6371.0
        >>> moon_mass =  CelestialDynamics.kg_to_solar_mass(7.34767309e22)
        >>> moon_distance = CelestialDynamics.km_to_AU(384400)
        >>> moon_velocity = CelestialDynamics.km_s_to_AU_yr(1.022)
        >>> moon_radius = 1737.0
        >>> earth = CelestialDynamics.celestial_object(earth_mass, [0,0,0], [0, 0, 0], colour='skyblue', radius=earth_radius)
        >>> moon = CelestialDynamics.celestial_object(moon_mass, [moon_distance, 0, 0], [0, moon_mass * moon_velocity, 0], colour='white', radius=moon_radius)
        '''
        self.mass = mass
        self.position = np.asarray(position, dtype=float)
        self.dot_position = np.array([0.0,0.0,0.0])
        self.momentum = np.asarray(momentum, dtype=float)
        self.dot_momentum = np.array([0.0,0.0,0.0])
        self.radius = radius
        self.colour = colour
        self.temp = temp
        # checks to make sure everything is the correct input type
        if self.position.shape != (3,) or self.momentum.shape != (3,):
            raise ValueError('Position and Momentum must be 1D arrays be of length 3.')

    def get_mass(self):
        return self.mass
    
    def get_pos(self):
        return self.position
    
    def get_vel(self):
        return self.dot_position
    
    def get_momentum(self):
        return self.momentum
    
    def get_force(self):
        return self.dot_momentum
    
    def get_colour(self):
        return self.colour
    
    def get_temp(self):
        return self.temp
    
    def get_radius(self):
        return self.radius
    
    def get_current_state(self):
        return self.mass, self.position, self.dot_position, self.momentum, self.dot_momentum, self.colour, self.temp
    
class stellar_system:
    def __init__(self, system: list):
        '''
        Bruh!

        Parameters:
            system (list): Collection of ``CelestialDynamics.celestial_object`` instances.

        Examples:
        >>> n = 10
        >>> stars = []

        >>> for i in range(n):
        >>>     mass_i = np.random.uniform(1, 10)
        >>>     pos_i = np.random.uniform(-3, 3, 3)
        >>>     momenta_i = mass_i * np.random.uniform(-10, 10, 3)
        >>>     stars.append(CelestialDynamics.celestial_object(mass_i, pos_i, momenta_i, mass_i)) 

        >>> system = CelestialDynamics.stellar_system(stars)
        '''
        self.system = system

    def get_stars(self):
        return self.system

    def get_state(self):
        positions = []
        momenta = []
        for star in self.system:
            q = star.get_pos()
            p = star.get_momentum()
            for i in range(3):
                positions.append(q[i])
                momenta.append(p[i])
        return *positions, *momenta

def integrate(system, t_eval, method='DOP853'):
    '''
    Bruh!

    Parameters:
        system (CelesticalDynamics.stellar_system): Collection of cellestial objects.

        t_eval (array_like): 1-D array as generated by ``np.linspace`` representing the period to integrate over.

        method (str, default = 'DOP853'): integration method to be used by ``solve_ivp``. See ``scipy.integrate``.

    Returns:
        output (scipy.integrate.solve_ivp): numerical solution to the system, and/or message from scipy integration. See ``scipy.integrate.solve_ivp``.
    Examples:
    >>> earth_mass = CelestialDynamics.kg_to_solar_mass(5.972e24)
    >>> earth_radius = 6371.0
    >>> moon_mass =  CelestialDynamics.kg_to_solar_mass(7.34767309e22)
    >>> moon_distance = CelestialDynamics.km_to_AU(384400)
    >>> moon_velocity = CelestialDynamics.km_s_to_AU_yr(1.022)
    >>> moon_radius = 1737.0

    >>> earth = CelestialDynamics.celestial_object(earth_mass, [0,0,0], [0, 0, 0], colour='skyblue', radius=earth_radius)
    >>> moon = CelestialDynamics.celestial_object(moon_mass, [moon_distance, 0, 0], [0, moon_mass * moon_velocity, 0], colour='white', radius=moon_radius)
    >>> system = CelestialDynamics.stellar_system([earth, moon])

    >>> t_eval = np.linspace(0, 1, 1001)
    >>> output = CelestialDynamics.integrate(system, t_eval, method='DOP853')
    >>> output
    message: The solver successfully reached the end of the integration interval.
    success: True
     status: 0
            t: [ 0.000e+00  1.000e-03 ...  9.990e-01  1.000e+00]
            y: [[ 0.000e+00  1.104e-07 ...  5.874e-05  5.749e-05]
                [ 0.000e+00  3.086e-09 ...  2.629e-03  2.634e-03]
                ...
                [ 7.959e-09  7.931e-09 ... -7.093e-09 -6.764e-09]
                [ 0.000e+00  0.000e+00 ...  0.000e+00  0.000e+00]]
         sol: None
    t_events: None
    y_events: None
        nfev: 647
        njev: 0
         nlu: 0
    '''
    def Hamiltonian_Vec_Field(t, PS):
        stars = system.get_stars()
        n = len(PS)//2
        dq = np.zeros([n])
        dp = np.zeros([n])
        n = n // 3
        for i in range(n):
            star1 = stars[i]
            for j in range(3):
                dq[3*i+j] = PS[3*(n+i)+j] / stars[i].get_mass() # dq/dt = DH/Dp = p/m
                pot = 0
                for k in range(len(stars)):
                    dist = 0
                    pot_t = 0
                    dist_i = 0
                    for l in range(3):
                        dist += (PS[3*i+l] - PS[3*k+l])**2 # compute the distance between stars
                    dist_i = (PS[3*i+j] - PS[3*k+j]) # compute distance along an axis
                    star2 = stars[k] #
                    if star2 != star1:
                        pot_t = star2.get_mass() * dist_i / dist**1.5 # calc the derivative of the potential wrt coordinate
                    pot += G * star1.get_mass() * pot_t # calc derivative of graviational potential energy 
                dp[3*i+j] = -pot # dp/dt = - DH/Dq = - DV/Dq
        return *dq, *dp
    t_span = [t_eval[0], t_eval[-1]]
    y0 = system.get_state()
    output = integ.solve_ivp(Hamiltonian_Vec_Field, t_span=t_span, y0=y0, t_eval=t_eval, method=method)
    return output

def animate_system_2D(system, t_eval, file_name, method='DOP853'):
    '''
    Bruh!

    Parameters:
        system (CelesticalDynamics.stellar_system): Collection of cellestial objects.

        t_eval (array_like): 1-D array as generated by ``np.linspace`` representing the period to integrate over.

        file_name (str): Name of the exported file. Must include file type (.mp4, .gif, ...).

        method (str, default = 'DOP853'): integration method to be used by ``solve_ivp``. See scipy.integrate.

    Examples:
    >>> earth_mass = CelestialDynamics.kg_to_solar_mass(5.972e24)
    >>> earth_radius = 6371.0
    >>> moon_mass =  CelestialDynamics.kg_to_solar_mass(7.34767309e22)
    >>> moon_distance = CelestialDynamics.km_to_AU(384400)
    >>> moon_velocity = CelestialDynamics.km_s_to_AU_yr(1.022)
    >>> moon_radius = 1737.0

    >>> earth = CelestialDynamics.celestial_object(earth_mass, [0,0,0], [0, 0, 0], colour='skyblue', radius=earth_radius)
    >>> moon = CelestialDynamics.celestial_object(moon_mass, [moon_distance, 0, 0], [0, moon_mass * moon_velocity, 0], colour='white', radius=moon_radius)
    >>> system = CelestialDynamics.stellar_system([earth, moon])

    >>> t_eval = np.linspace(0, 1, 1001)
    >>> CelestialDynamics.animate_system_2D(system, t_eval, 'earth_moon.mp4')
    '''
    stars = system.get_stars()
    output = integrate(system ,t_eval, method=method)
    plt.rcParams['figure.facecolor'] = 'black'
    frames = len(t_eval)
    inter = 1000/frames
    fig = plt.figure()
    ax = plt.axes(aspect='equal')
    ax.set_axis_off()
    min_axis = 0
    max_axis = 0
    max_radius = 0
    for k in range(len(stars)):
        for l in range(3):
            if np.min(output.y[3*k + l]) < min_axis:
                min_axis = np.min(output.y[3*k + l])
            if np.max(output.y[3*k + l]) > max_axis:
                max_axis = np.max(output.y[3*k + l])
        if stars[k].get_radius() > max_radius:
                max_radius = stars[k].get_radius()
        if np.abs(max_axis) > np.abs(min_axis):
            min_axis = -max_axis
        else:
            max_axis = -min_axis
    def animate(i):
        ax.clear()
        ax.set_axis_off()
        ax.set_xlim(1.5*min_axis, 1.5*max_axis)
        ax.set_ylim(1.5*min_axis, 1.5*max_axis)
        for j in range(len(stars)):
            plt.scatter(output.y[3*j][i], output.y[3*j+1][i], color=stars[j].get_colour(), s = 50*stars[j].get_radius()/max_radius)
            plt.plot(output.y[3*j][0:i], output.y[3*j+1][0:i], color=stars[j].get_colour(), alpha=0.5)
    ani = animation.FuncAnimation(fig, animate, frames=frames, interval=inter, repeat=True)
    ffmpeg_writer = animation.FFMpegWriter(fps=60)
    print('Rendering Animation:')
    ani.save(file_name, writer=ffmpeg_writer)
    print('Animation Rendered Successfully.')
