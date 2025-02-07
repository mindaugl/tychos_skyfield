"""
Test script showing how to set up multiple object(planet) system.
The test agrees with the Tychosium - both planet locations and quaternions were tested.
"""
from scipy.spatial.transform import Rotation as R
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
        radius_vec
        children

    Methods
    -------
        __init__
        move_planet
        add_child

    Notes
    -----
    move_planet() method needs to be called for parent first and only then for child.
    """
    def __init__(self, orbit_radius=100.0, orbit_center = OrbitCenter(),
                 orbit_tilt = OrbitTilt(), start_pos=20.0):

        self.start_pos = start_pos
        self.rotation = (R.from_euler('x', orbit_tilt.x, degrees = True) *
                          R.from_euler('z', orbit_tilt.z, degrees = True))
        self.location = np.array([0.0, 0.0, 0.0])
        self.center = (np.array([orbit_center.x, orbit_center.y, orbit_center.z]).
                       astype(np.float64))
        self.radius_vec = np.array([orbit_radius, 0.0, 0.0])
        self.children = []

    def move_planet(self, speed, pos):
        """
        Moves planet by specified pos and speed
        :param speed: float
            Moving speed multiplier.
        :param pos: float
            Position in degrees to rotate around y-axis
        :return:
        """
        self.rotation = self.rotation * R.from_euler('y', speed * pos - self.start_pos,
                                                     degrees = True)
        radius_rotated = self.rotation.apply(self.radius_vec)
        self.location = self.center + radius_rotated
        for child in self.children:
            child.rotation = self.rotation * child.rotation
            child.center = self.center + self.rotation.apply(self.radius_vec + child.center)

    def add_child(self, child_obj):
        """
        Add child to the planet.
        Order of move_planet() matters, need to move parent first
        :param child_obj: PlanetObj
            Child object to be added.
        :return: none
        """
        self.children += [child_obj]

### Check for 2 planet system, parameters arbitrary for full parameter scope coverage,
# both quaternions and positions agree with tychosium

sun_temp_def = PlanetObj(35, OrbitCenter(10, 20, -10),
                    OrbitTilt(30, 60), 180)
sun_temp = PlanetObj(100, OrbitCenter(25, 35, -10),
                OrbitTilt(15, 30), 15)
sun_temp_def.add_child(sun_temp)
sun_temp_def.move_planet(1, 0)
print("sun_temp before move:")
print("location:", sun_temp.location)
print("rotation vector:", sun_temp.rotation.as_rotvec())
print("quaternion:", sun_temp.rotation.as_quat())
sun_temp.move_planet(1, 0)

print("sun_temp_def:")
print("location:", sun_temp_def.location)
print("rotation vector:", sun_temp_def.rotation.as_rotvec())
print("quaternion:", sun_temp_def.rotation.as_quat())

print("sun_temp:")
print("location:", sun_temp.location)
print("rotation vector:", sun_temp.rotation.as_rotvec())
print("quaternion:", sun_temp.rotation.as_quat())

### Check mercury planet system, parameters are actual from tychosium,
# both quaternions and positions agree with tychosium

earth = PlanetObj( 37.8453, OrbitCenter(0, 0, 0),
                      OrbitTilt(0, 0), 0)
mer_def_a = PlanetObj(100, OrbitCenter(-6.9, -3.2, 0),
                      OrbitTilt(0, 0), 0)
mer_def_b = PlanetObj(0, OrbitCenter(0, 0, 0),
                      OrbitTilt(-1.3, 0.5), 33)
mer = PlanetObj(38.710225, OrbitCenter(0.6, 3, -0.1),
                OrbitTilt(3, 0.5), -180.8)

earth.add_child(mer_def_a)
mer_def_a.add_child(mer_def_b)
mer_def_b.add_child(mer)

position = 0 # the new position in degrees

earth.move_planet(1, 90 + position)
mer_def_a.move_planet(1, position)
mer_def_b.move_planet(1, position)
mer.move_planet(1, position)
print("")

print("earth:")
print("location:", earth.location)
print("rotation vector:", earth.rotation.as_rotvec())
print("quaternion:", earth.rotation.as_quat())

print("mer_def_a:")
print("location:", mer_def_a.location)
print("rotation vector:", mer_def_a.rotation.as_rotvec())
print("quaternion:", mer_def_a.rotation.as_quat())

print("mer_def_b:")
print("location:", mer_def_b.location)
print("rotation vector:", mer_def_b.rotation.as_rotvec())
print("quaternion:", mer_def_b.rotation.as_quat())

print("mer:")
print("location:", mer.location)
print("rotation vector:", mer.rotation.as_rotvec())
print("quaternion:", mer.rotation.as_quat())

