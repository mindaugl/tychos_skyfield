"""
Test script showing how to set up 2 object(planet) system - parent and child.
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

sun_def = PlanetObj(35, OrbitCenter(10, 20, -10),
                    OrbitTilt(30, 60), 180)
sun = PlanetObj(100, OrbitCenter(25, 35, -10),
                OrbitTilt(15, 30), 15)
sun_def.add_child(sun)
sun_def.move_planet(1, 0)
print("sun before move:")
print("location:", sun.location)
print("rotation vector:", sun.rotation.as_rotvec())
print("quaternion:", sun.rotation.as_quat())
sun.move_planet(1, 0)

print("sun_def:")
print("location:", sun_def.location)
print("rotation vector:", sun_def.rotation.as_rotvec())
print("quaternion:", sun_def.rotation.as_quat())

print("sun:")
print("location:", sun.location)
print("rotation vector:", sun.rotation.as_rotvec())
print("quaternion:", sun.rotation.as_quat())
