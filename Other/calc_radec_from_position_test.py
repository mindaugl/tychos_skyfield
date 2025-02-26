"""
Test module for calculting radec from the position vector
"""
from skyfield.api import load
from skyfield.positionlib import Astrometric
import numpy as np

def radec_from_position(observer, target_pos, target_vel, t):
    """
    Calculate ra, dec, distance using target position vector as seen from observer.
    No time of light travel or other adjustments applied.
    :param observer: observer object used as reference frame (VectorSum)
    :param target_pos: np array for xyz coordinates for target position (au)
    :param target_vel: np array for target velocities (au/d) in xyz coordinates,
        currently not relevant
    :param t: time when the positions are taken, instantaneously (Timescale)
    :return: ra, dec, distance from observer to target
    """
    baryocentric = observer.at(t)
    obs_position = baryocentric.position.au
    obs_velocity = baryocentric.velocity.au_per_d
    astrometric = Astrometric(target_pos - obs_position, target_vel - obs_velocity, time)
    return astrometric.radec()

planets = load('de421.bsp')
obs = planets["earth"]
ts = load.timescale()
time = ts.utc(2000,6,21,12,0,0)

### Test arbitrary position:
t_position = np.array([1, 1, 1])# au
t_velocity = np.array([0, 0, 0]) # au/d
ra, dec, distance = radec_from_position(obs, t_position, t_velocity, time)
print("Test arb position: ", ra, dec, distance)

### Test planet as calculated from ephemeris vs position:
name_planet = "NEPTUNE_BARYCENTER"
planet = planets[name_planet]
baryo = obs.at(time)
astrom = baryo.observe(planet)
ra, dec, distance = astrom.radec()
print("\nTest", name_planet, "as calc by ephemeris:  ", ra, dec, distance)
baryo = planet.at(time)
t_position = baryo.position.au
t_velocity = baryo.velocity.au_per_d
ra, dec, distance = radec_from_position(obs, t_position, t_velocity, time)
print("Test", name_planet, "as calc from position: ", ra, dec, distance)
