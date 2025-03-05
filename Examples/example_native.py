"""
Example script showing the usage of Tychos model in Python3.
Library dependencies are: numpy and scipy.
"""

import datetime as dt
from tychos_skyfield import baselib as T


def get_julian_day(date_time):
    """
    Helper function to obtain Julian Day for a given datetime
    :param date_time: datetime
        UCT time
    :return: float
    """
    ref = dt.datetime(2000,6,21,12,0,0)
    julian_day = (date_time - ref).total_seconds()/24/3600 + 2451717.0
    return julian_day

# Create solar system:
system = T.TychosSystem()

# move system to particular time:
time = dt.datetime(2020, 6, 21, 12, 0, 0)
system.move_system(get_julian_day(time))

# Check list of observable objects:
print("Observable objects:", system.get_observable_objects())

# Check location of object, given in centi-AU (earth-sun distance is ~100):
print("Location:", system['Jupiter'].location)

# Get RA/Dec/Distance as calculated in Tychosium (in the frame associated with Earth's orientation):
ra, dec, dist = system['Jupiter'].radec_direct(system['Earth'], system['Polar_axis'], 'date')
print("RA/Dec/Distance in moving frame:   ", ra, ",", dec, ",", dist)

# Get RA/Dec/Distance for the J2000 epoch (in the frame of Earth being on the 2000/01/01):
ra, dec, dist = system['Jupiter'].radec_direct(system['Earth'], None, 'j2000')
print("RA/Dec/Distance in the j2000 frame:", ra, ",", dec, ",", dist)
