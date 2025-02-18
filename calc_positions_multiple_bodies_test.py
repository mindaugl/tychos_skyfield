"""
Test script showing how to set up multiple object(planet) system.
Can specify the datetime or degrees for how to move the planets.
Calculates RA and DEC of planet directly.
The test agrees with the Tychosium.

"""
from scipy.spatial.transform import Rotation as R
from skyfield.api import load
import numpy as np

class OrbitCenter:
    """
    Data class to keep orbit center coordinates
    """
    def __init__(self, orbit_center_a = 0.0, orbit_center_b = 0.0, orbit_center_c = 0.0):
        self.x = orbit_center_a
        self.y = orbit_center_c
        self.z = orbit_center_b

class OrbitTilt:
    """
    Data class to keep orbit tilt values
    """
    def __init__(self, orbit_tilt_a=30.0, orbit_tilt_b=60.0):
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

    Notes
    -----
    move_planet() method needs to be called for parent first and only then for child.
    speed = 1/period(years) * 2pi, represents rotation radians per year
    """

    def __init__(self, orbit_radius = 100.0, orbit_center = OrbitCenter(),
                 orbit_tilt = OrbitTilt(), start_pos = 20.0, speed = 0.0):

        self.start_pos = start_pos
        self.rotation = (R.from_euler('x', orbit_tilt.x, degrees = True) *
                          R.from_euler('z', orbit_tilt.z, degrees = True))
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
        pos = (time.tt - 2451717.0)/365.2425 * 360
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

    def move_planet_basic(self, pos):
        """
        Moves planet by specified pos, assuming self.speed = 0 and self.start_pos = 0.
        Can call this function multiple times - it does not modify children.
        :param pos: float
            Position in degrees to rotate around y-axis
        :return: none
        """
        self.rotation = self.rotation * R.from_euler('y', pos, degrees = True)
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

    def radec_direct(self):
        """
        Calculate RA and DEC for the current location of the planet. It uses direct coordinate
        transformation to the earth rotational axis/plane to calculate RA and DEC.
        :return: tuple - (ra, dec)
            ra is calculated in hours (and fractions of hours)
            dec is calculated in degrees (and fractions of degrees)
        """
        rot = R.from_euler('zxy', [-23.439062, 0.26, 90], degrees = True)
        unit_prime = rot.apply(np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
        loc_prime = np.dot(unit_prime, self.location)
        dec = 90 - (np.arccos(loc_prime[1] / np.sqrt(np.dot(loc_prime, loc_prime))) * 180 / np.pi)
        ra = (np.sign(loc_prime[0]) *
              np.arccos(loc_prime[2] / np.sqrt(loc_prime[0] ** 2 + loc_prime[2] ** 2)))
        ra = ra * 12 / np.pi
        return (ra, dec)

### Check for 2 planet system, parameters arbitrary for full parameter scope coverage,
# both quaternions and positions agree with tychosium

sun_temp_def = PlanetObj(35, OrbitCenter(10, 20, -10),
                    OrbitTilt(30, 60), 180, 0)
sun_temp = PlanetObj(100, OrbitCenter(25, 35, -10),
                OrbitTilt(15, 30), 15, 0)
sun_temp_def.add_child(sun_temp)
sun_temp_def.move_planet(0)
print("sun_temp before move:")
print("location:", sun_temp.location)
print("rotation vector:", sun_temp.rotation.as_rotvec())
print("quaternion:", sun_temp.rotation.as_quat())
sun_temp.move_planet(0)

print("sun_temp_def:")
print("location:", sun_temp_def.location)
print("rotation vector:", sun_temp_def.rotation.as_rotvec())
print("quaternion:", sun_temp_def.rotation.as_quat())

print("sun_temp:")
print("location:", sun_temp.location)
print("rotation vector:", sun_temp.rotation.as_rotvec())
print("quaternion:", sun_temp.rotation.as_quat())

### Check mercury planet system, parameters are actual from tychosium,
# use position of degrees,
# both quaternions and positions agree with tychosium

earth = PlanetObj( 37.8453, OrbitCenter(0, 0, 0),
                      OrbitTilt(0, 0), 0, -0.0002479160869310127)
mer_def_a = PlanetObj(100, OrbitCenter(-6.9, -3.2, 0),
                      OrbitTilt(0, 0), 0, 2 * np.pi)
mer_def_b = PlanetObj(0, OrbitCenter(0, 0, 0),
                      OrbitTilt(-1.3, 0.5), 33, -2 * np.pi)
mer = PlanetObj(38.710225, OrbitCenter(0.6, 3, -0.1),
                OrbitTilt(3, 0.5), -180.8, 26.08763045)

earth.add_child(mer_def_a)
mer_def_a.add_child(mer_def_b)
mer_def_b.add_child(mer)

position = 0.25 * 360 # the new position in degrees. For speed = 2pi, position repeats at pos = 360

earth.move_planet_basic(90)
earth.move_planet(position)
mer_def_a.move_planet(position)
mer_def_b.move_planet(position)
mer.move_planet(position)
print("")

print("earth:")
print("location:", earth.location)
print("rotation vector:", earth.rotation.as_rotvec())
print("quaternion:", earth.rotation.as_quat())

print("mercury_def_a:")
print("location:", mer_def_a.location)
print("rotation vector:", mer_def_a.rotation.as_rotvec())
print("quaternion:", mer_def_a.rotation.as_quat())

print("mercury_def_b:")
print("location:", mer_def_b.location)
print("rotation vector:", mer_def_b.rotation.as_rotvec())
print("quaternion:", mer_def_b.rotation.as_quat())

print("mercury:")
print("location:", mer.location)
print("rotation vector:", mer.rotation.as_rotvec())
print("quaternion:", mer.rotation.as_quat())

### check using datetime

earth = PlanetObj( 37.8453, OrbitCenter(0, 0, 0),
                      OrbitTilt(0, 0), 0, -0.0002479160869310127)
mer_def_a = PlanetObj(100, OrbitCenter(-6.9, -3.2, 0),
                      OrbitTilt(0, 0), 0, 2 * np.pi)
mer_def_b = PlanetObj(0, OrbitCenter(0, 0, 0),
                      OrbitTilt(-1.3, 0.5), 33, -2 * np.pi)
mer = PlanetObj(38.710225, OrbitCenter(0.6, 3, -0.1),
                OrbitTilt(3, 0.5), -180.8, 26.08763045)

earth.add_child(mer_def_a)
mer_def_a.add_child(mer_def_b)
mer_def_b.add_child(mer)

ts = load.timescale()
t = ts.tt(2024,9,21,16,0,30)

earth.move_planet_basic(90)
earth.move_planet_tt(t)
mer_def_a.move_planet_tt(t)
mer_def_b.move_planet_tt(t)
mer.move_planet_tt(t)

print("earth:")
print("location:", earth.location)
print("rotation vector:", earth.rotation.as_rotvec())
print("quaternion:", earth.rotation.as_quat())

print("mercury_def_a:")
print("location:", mer_def_a.location)
print("rotation vector:", mer_def_a.rotation.as_rotvec())
print("quaternion:", mer_def_a.rotation.as_quat())

print("mercury_def_b:")
print("location:", mer_def_b.location)
print("rotation vector:", mer_def_b.rotation.as_rotvec())
print("quaternion:", mer_def_b.rotation.as_quat())

print("mercury:")
print("location:", mer.location)
print("rotation vector:", mer.rotation.as_rotvec())
print("quaternion:", mer.rotation.as_quat())

### Check RA DEC calculation, tested to agree with Tychosium

earth = PlanetObj( 0.0, OrbitCenter(0, 0, 0),
                      OrbitTilt(0, 0), 0, -0.0002479160869310127)
# note that Earth is set to be immovable and located at (0, 0, 0)
mer_def_a = PlanetObj(100, OrbitCenter(-6.9, -3.2, 0),
                      OrbitTilt(0, 0), 0, 2 * np.pi)
mer_def_b = PlanetObj(0, OrbitCenter(0, 0, 0),
                      OrbitTilt(-1.3, 0.5), 33, -2 * np.pi)
mer = PlanetObj(38.710225, OrbitCenter(0.6, 3, -0.1),
                OrbitTilt(3, 0.5), -180.8, 26.08763045)

earth.add_child(mer_def_a)
mer_def_a.add_child(mer_def_b)
mer_def_b.add_child(mer)

position = 0.0

earth.move_planet_basic(90)
earth.move_planet(position)
mer_def_a.move_planet(position)
mer_def_b.move_planet(position)
mer.move_planet(position)

(ra_, dec_) = mer.radec_direct()

print("\nmercury declination (degrees):", dec_, ", RA (hours):", ra_)

