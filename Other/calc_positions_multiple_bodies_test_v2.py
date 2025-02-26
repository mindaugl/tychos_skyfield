"""
Tychosium model implementation in Python.
Can specify the datetime or degrees for how to move the planets.
Calculates RA and DEC of planet directly and using Skyfield Astrometric position.
Compares RA, DEC with ones obtained in Skyfield for main objects.
Uses updated Tychosium code as reference.
"""
from scipy.spatial.transform import Rotation as R
from skyfield.api import load, Angle
from skyfield.positionlib import Astrometric
import numpy as np


class OrbitCenter:
    """
    Data class to keep orbit center coordinates
    """

    def __init__(self, orbit_center_a=0.0, orbit_center_b=0.0, orbit_center_c=0.0):
        self.x = orbit_center_a
        self.y = orbit_center_c
        self.z = orbit_center_b


class OrbitTilt:
    """
    Data class to keep orbit tilt values
    """

    def __init__(self, orbit_tilt_a=0.0, orbit_tilt_b=0.0):
        self.x = orbit_tilt_a
        self.z = orbit_tilt_b


class PlanetObj:
    """
    Class for planet object definition.
    It initializes the starting parameters, calculates the rotations, and locations.

    Attributes
    ----------
        start_pos
        rotation
        location
        center
        speed
        radius_vec
        children

    Methods
    -------
        __init__
        move_planet_tt
        move_planet
        move_planet_basic
        add_child
        radec_direct
        radec

    Notes
    -----
    move_planet() method needs to be called for parent first and only then for child.
    speed = 1/period(years) * 2pi, represents rotation radians per year
    """

    def __init__(self, orbit_radius=100.0, orbit_center=OrbitCenter(),
                 orbit_tilt=OrbitTilt(), start_pos=20.0, speed=0.0):

        self.start_pos = start_pos
        self.rotation = (R.from_euler('x', orbit_tilt.x, degrees=True) *
                         R.from_euler('z', orbit_tilt.z, degrees=True))
        self.location = np.array([0.0, 0.0, 0.0])
        self.center = (np.array([orbit_center.x, orbit_center.y, orbit_center.z]).
                       astype(np.float64))
        self.radius_vec = np.array([orbit_radius, 0.0, 0.0])
        self.speed = speed / (2 * np.pi)
        self.children = []

    def move_planet_tt(self, time):
        """
        Moves planet to specified terrestrial datetime time.
        NOTE: only can use function once, as everytime it modifies children values.
        :param time: Time
            Skyfield Time terrestrial time value to which to move the planet
        :return: none
        """
        pos = (time.tt - 2451717.0) / 365.2425 * 360
        # 2451717 is reference Julian Date tt for date 2000-6-21 12:00:00
        self.move_planet(pos)

    def move_planet(self, pos):
        """
        Moves planet by specified pos.
        NOTE: only can use function once, as everytime it modifies children values.
        :param pos: float
            Position in degrees to rotate around y-axis
        :return: none
        """
        self.move_planet_basic(self.speed * pos - self.start_pos)
        for child in self.children:
            child.rotation = self.rotation * child.rotation
            child.center = self.center + self.rotation.apply(self.radius_vec + child.center)

    def move_planet_basic(self, pos, directions='y'):
        """
        Moves planet by specified pos, assuming self.speed = 0 and self.start_pos = 0.
        Can call this function multiple times - it does not modify children.
        :param pos: float or List[float]
            Position(s) in degrees to rotate around 'directions'
        :param directions: [optional] string
            The direction or multiple directions with respect which to move
        :return: none
        """
        self.rotation = self.rotation * R.from_euler(directions, pos, degrees=True)
        radius_rotated = self.rotation.apply(self.radius_vec)
        self.location = self.center + radius_rotated

    def add_child(self, child_obj):
        """
        Add child to the planet.
        Order of move_planet() matters, need to move parent first
        :param child_obj: PlanetObj
            Child object to be added.
        :return: none
        """
        self.children += [child_obj]

    def radec_direct(self, ref_obj, polar_obj, epoch='j2000'):
        """
        Calculate RA and DEC for the current location of the planet. It uses projects planet
        location to the appropriate ref frame for the epoch
        :param ref_obj: PlanetObj
            reference object with respect to which calculate RA and DEC, typically earth
        :param polar_obj: PlanetObj
            reference object that contains transformation for polar axis frame which
            is used to calculate RA, DEC
            Only required for the epoch = 'date'
        :param epoch: Optional[String]: 'j2000'(default) or 'date'
            epoch specifies which 'time' is used for ra/dec calculation. 'j2000' corresponds
            to J2000(or ICRF), while 'date' is frame associated with current time
        :return: tuple[Angle, Angle, Float] - (ra, dec, dist)
            ra is calculated in hours (and fractions of hours)
            dec is calculated in degrees (and fractions of degrees)
            dist is the distance to the planet from earth in AU
        """
        if epoch == 'j2000':
            rot = R.from_euler('zxy', [-23.439062, 0.26, 90], degrees=True)
        elif epoch == 'date':
            rot = polar_obj.rotation
        else:
            raise AttributeError("Unknown epoch provided: " + epoch +
                            ". Only epochs 'j2000' and 'date' are supported." )

        unit_prime = rot.apply(np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
        loc_prime = np.dot(unit_prime, self.location - ref_obj.location)
        dec = np.pi / 2 - np.arccos(loc_prime[1] / np.sqrt(np.dot(loc_prime, loc_prime)))
        dec = Angle(radians=dec, preference='degrees', signed=True)
        ra = (np.sign(loc_prime[0]) *
              np.arccos(loc_prime[2] / np.sqrt(loc_prime[0] ** 2 + loc_prime[2] ** 2)))
        if ra < 0:
            ra += 2 * np.pi
        ra = Angle(radians=ra, preference='hours')
        dist = np.linalg.norm(loc_prime) / 100
        return ra, dec, dist

    def radec(self, time, ref_obj, polar_obj, epoch='j2000'):
        """
        Calculate ra, dec, distance using target position vector as seen from observer
        using Skyfield Astrometric position.
        No time of light travel adjustments applied (that are applied in Skyfield by default).
        Agrees with radec_direct() function.
        :param time: Time
            time when the positions and radec are calculated, instantaneous
        :param ref_obj: PlanetObj
            reference object with respect to which calculate RA and DEC, typically earth
        :param polar_obj: PlanetObj
            reference object that contains transformation for polar axis frame to
            which Tychos coordinates are transformed to pass to the Skyfield
            Only required for the epoch = 'date'
        :param epoch: Optional[String]: 'j2000'(default) or 'date'
            epoch specifies which 'time' is used for ra/dec calculation. 'j2000' corresponds
            to J2000(or ICRF), while 'date' is frame associated with current time
        :return: Tuple[Angle, Angle, Distance]
            ra, dec, distance from observer to target
        """
        if epoch == 'j2000':
            r1 = R.from_euler('ZXY', [-23.439062, 0.26, 90], degrees=True) # intrinsic rotation
            r2 = R.from_euler('yx', [-90, 90], degrees=True) # extrinsic to align axis with Skyfield
        elif epoch == 'date':
            r1 = polar_obj.rotation.inv()
            r2 = R.from_euler('zy', [90, 90], degrees=True) # extrinsic to align axis with Skyfield
        else:
            raise AttributeError("Unknown epoch provided: " + epoch +
                            ". Only epochs 'j2000' and 'date' are supported." )

        loc_new = (r2 * r1).apply(self.location - ref_obj.location) / 100
        astrometric = Astrometric(loc_new, np.array([0, 0, 0]), time)
        return astrometric.radec()

##### Compare to actual Skyfield positions, neglecting light travel time corrections ####

ts = load.timescale()
t = ts.tt(1000, 6, 21, 12, 0, 0)

# planets = load('de421.bsp')
planets = load('de422.bsp')
obs = planets["earth"]
barycentric = obs.at(t)
print("\nt:", t.tt_strftime())

earth = PlanetObj(37.8453, OrbitCenter(0, 0, 0),
                  OrbitTilt(0, 0), 0, -0.0002479160869310127)
polar_axis = PlanetObj(0, OrbitCenter(0, 0, 0),
                       OrbitTilt(0, 0), 0, 0.0)

sun_def = PlanetObj(0.0, OrbitCenter(1.4, -0.6, 0.0),
                    OrbitTilt(0.1, 0.0), 0.0, 0.0)
sun = PlanetObj(100.0, OrbitCenter(1.2, -0.1, 0.0),
                OrbitTilt(0.1, 0.0), 0.0, 2 * np.pi)

mercury_def_a = PlanetObj(100, OrbitCenter(-6.9, -3.2, 0),
                          OrbitTilt(0, 0), 0, 2 * np.pi)
mercury_def_b = PlanetObj(0, OrbitCenter(0, 0, 0),
                          OrbitTilt(-1.3, 0.5), 33, -2 * np.pi)
mercury = PlanetObj(38.710225, OrbitCenter(0.6, 3, -0.1),
                    OrbitTilt(3, 0.5), -180.8, 26.08763045)

m_factor = 39.2078
moon_def_a = PlanetObj(0.0279352315075 / m_factor,
                       OrbitCenter(0 / m_factor, 0 / m_factor, 0 / m_factor),
                       OrbitTilt(-0.2, 0.5), 226.4, 0.71015440177343)
moon_def_b = PlanetObj(0 / m_factor, OrbitCenter(-0.38 / m_factor, 0.22 / m_factor, 0 / m_factor),
                       OrbitTilt(2.3, 2.6), -1.8, 0.0)
moon = PlanetObj(10 / m_factor, OrbitCenter(0.8 / m_factor, -0.81 / m_factor, -0.07 / m_factor),
                 OrbitTilt(-1.8, -2.6), 261.2, 83.28521)

venus_def_a = PlanetObj(100, OrbitCenter(0.5, 0.5, 0),
                        OrbitTilt(0, 0), 0, 2 * np.pi)
venus_def_b = PlanetObj(0, OrbitCenter(0, 0.65, 0),
                        OrbitTilt(0, 0), 16.6, -2 * np.pi)
venus = PlanetObj(72.327789, OrbitCenter(0.6, -0.9, 0),
                  OrbitTilt(3.2, -0.05), -23.6, 10.21331385)

mars_def_e = PlanetObj(100, OrbitCenter(10.1, -20.7, 0),
                       OrbitTilt(0, 0), 0, 2 * np.pi)
mars_def_s = PlanetObj(7.44385, OrbitCenter(0, 0, 0),
                       OrbitTilt(0, 0), -115, 0.3974599)
mars = PlanetObj(152.677, OrbitCenter(0, 0, 0),
                 OrbitTilt(-0.2, -1.7), 119.3, -3.33985)

phobos = PlanetObj(5, OrbitCenter(0, 0, 0),
                   OrbitTilt(0, 0), 122, 6986.5)
deimos = PlanetObj(10, OrbitCenter(0, 0, 0),
                   OrbitTilt(0, 0), 0, 1802.0)

jupiter_def = PlanetObj(0.0, OrbitCenter(0.0, 0.0, 0.0),
                        OrbitTilt(0.0, 0.0), 75.4, -2 * np.pi)
jupiter = PlanetObj(520.4, OrbitCenter(-49.0, 3.0, -1.0), OrbitTilt(0.0, -1.2),
                    -34.0, 0.52994136)

saturn_def = PlanetObj(20, OrbitCenter(11, 0, 0),
                       OrbitTilt(0, 0), 518, -2 * np.pi)
saturn = PlanetObj(958.2, OrbitCenter(69, 40, 0),
                   OrbitTilt(-2.5, 0), -123.8, 0.21351984)

uranus_def = PlanetObj(20, OrbitCenter(0, 0, 0),
                       OrbitTilt(0, 0), 123, -2 * np.pi)
uranus = PlanetObj(1920.13568, OrbitCenter(150, -65, 0),
                   OrbitTilt(-0.2, -0.7), 371.8, 0.07500314)

neptune_def = PlanetObj(20, OrbitCenter(0, 0, 0),
                        OrbitTilt(0, 0), 175.2, -2 * np.pi)
neptune = PlanetObj(3004.72, OrbitCenter(0, 20, 0),
                    OrbitTilt(-1.6, 1.15), 329.3, 0.03837314)

halleys_def = PlanetObj(20, OrbitCenter(-5, 10, 11),
                        OrbitTilt(0, 0), 179, -2 * np.pi)
halleys = PlanetObj(1674.5, OrbitCenter(-1540, -233.5, -507),
                    OrbitTilt(6.4, 18.55), 76.33, -0.0830100973)

eros_def_a = PlanetObj(100, OrbitCenter(-40, 31.5, -0.5),
                       OrbitTilt(-7.3, 3.6), 0, 2 * np.pi)
eros_def_b = PlanetObj(0, OrbitCenter(-16, -4.5, 0),
                       OrbitTilt(0, 0), 0, -7.291563307179587)
eros = PlanetObj(145.79, OrbitCenter(5.2, -6, 0),
                 OrbitTilt(0, 0), 171.8, 4.57668492)

d_pl = {'earth': earth, 'polar_axis': polar_axis, 'sun_def': sun_def, 'sun': sun,
        'mercury_def_a': mercury_def_a, 'mercury_def_b': mercury_def_b, 'mercury': mercury,
        'moon_def_a': moon_def_a, 'moon_def_b': moon_def_b, 'moon': moon,
        'venus_def_a': venus_def_a, 'venus_def_b': venus_def_b, 'venus': venus,
        'mars_def_e': mars_def_e, 'mars_def_s': mars_def_s, 'mars': mars,
        'phobos': phobos, 'deimos': deimos, 'jupiter_def': jupiter_def, 'jupiter': jupiter,
        'saturn_def': saturn_def, 'saturn': saturn,
        'uranus_def': uranus_def, 'uranus': uranus,
        'neptune_def': neptune_def, 'neptune': neptune,
        'halleys_def': halleys_def, 'halleys': halleys,
        'eros_def_a': eros_def_a, 'eros_def_b': eros_def_b, 'eros': eros}

all_planets = ['earth', 'polar_axis', 'sun_def', 'sun', 'mercury_def_a', 'mercury_def_b', 'mercury',
               'moon_def_a', 'moon_def_b', 'moon', 'venus_def_a', 'venus_def_b', 'venus',
               'mars_def_e', 'mars_def_s', 'mars', 'phobos', 'deimos', 'jupiter_def', 'jupiter',
               'saturn_def', 'saturn', 'uranus_def', 'uranus', 'neptune_def', 'neptune',
               'halleys_def', 'halleys', 'eros_def_a', 'eros_def_b', 'eros']

print_planets = ['sun', 'mercury', 'moon', 'venus', 'mars', 'phobos', 'deimos',
                 'jupiter', 'saturn', 'uranus', 'neptune', 'halleys', 'eros']

skyfield_map = {'sun': 'SUN', 'moon': 'MOON', 'earth': 'EARTH', 'mercury': 'MERCURY',
                'venus': 'VENUS', 'mars': 'MARS', 'jupiter': 'JUPITER_BARYCENTER',
                'saturn': 'SATURN_BARYCENTER', 'uranus': 'URANUS_BARYCENTER',
                'neptune': 'NEPTUNE_BARYCENTER', 'pluto': 'PLUTO_BARYCENTER'}

earth.add_child(polar_axis)

earth.add_child(sun_def)
sun_def.add_child(sun)

earth.add_child(moon_def_a)
moon_def_a.add_child(moon_def_b)
moon_def_b.add_child((moon))

earth.add_child(mercury_def_a)
mercury_def_a.add_child(mercury_def_b)
mercury_def_b.add_child(mercury)

earth.add_child(venus_def_a)
venus_def_a.add_child(venus_def_b)
venus_def_b.add_child(venus)

earth.add_child(mars_def_e)
mars_def_e.add_child(mars_def_s)
mars_def_s.add_child(mars)

mars.add_child(phobos)
mars.add_child(deimos)

sun.add_child(jupiter_def)
jupiter_def.add_child(jupiter)

sun.add_child(saturn_def)
saturn_def.add_child(saturn)

sun.add_child(uranus_def)
uranus_def.add_child(uranus)

sun.add_child(neptune_def)
neptune_def.add_child(neptune)

sun.add_child(halleys_def)
halleys_def.add_child(halleys)

earth.add_child(eros_def_a)
eros_def_a.add_child(eros_def_b)
eros_def_b.add_child(eros)

polar_axis.move_planet_basic([-23.439062, 0.26], 'zx')
earth.move_planet_basic(90)

for p in all_planets:
    d_pl[p].move_planet_tt(t)

for p in print_planets:
    print("\n", p)
    print("Tyhos             :", d_pl[p].radec(t, earth, polar_axis, 'date'))
    print("Tyhos direct      :", d_pl[p].radec_direct(earth, polar_axis, 'date'))
    print("Tyhos j2000       :", d_pl[p].radec(t, earth, 'j2000'))
    print("Tyhos direct j2000:", d_pl[p].radec_direct(earth, 'j2000'))
    if p in skyfield_map:
        print("Skyfield epoch 'date':", barycentric.observe(planets[skyfield_map[p]]).radec('date'))
        print("Skyfield epoch ICRF  :", barycentric.observe(planets[skyfield_map[p]]).radec())
        print("Skyfield epoch J2000 :",
              barycentric.observe(planets[skyfield_map[p]]).radec(ts.J2000))
    print("location:", d_pl[p].location)
    print("quaternion:", d_pl[p].rotation.as_quat())
