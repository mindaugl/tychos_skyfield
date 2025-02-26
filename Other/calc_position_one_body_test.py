"""
Example of calculating the static position and quaterion of the object (planet)
    that has no deferents.
This was checked to agrees with the Tychosium simulation that uses three.js.
"""
from scipy.spatial.transform import Rotation as R
import numpy as np

### set up one object (Sun), parameter names taken from tychosium:
orbitRadius = 100 # 1 au
orbitCentera = 1.2 # centi-au
orbitCenterb = -0.1 # centi-au
orbitCenterc = 0 # centi-au
orbitTilta = 30 # degrees
orbitTiltb = 60 # degrees
startPos = 20 # starting angle, degrees
speed = 1 # degrees per time
pos = 0 # current angle, degrees

### apply tilt around x and z axis and then rotate around y axis:
r_orbit = R.from_euler('XZY', [orbitTilta, orbitTiltb, pos * speed - startPos],
                       degrees = True) # orbit rotation is intrinsic
planet = np.array([orbitRadius, 0, 0])
planet_2 = r_orbit.apply(planet)
planet_3 = planet_2 + np.array([orbitCentera, orbitCenterc, orbitCenterb])

print("Planet position at pos", pos, ":", planet_3)
print("Planet quaterion at pos", pos, ":", r_orbit.as_quat())
