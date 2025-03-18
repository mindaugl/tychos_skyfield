"""
Check radec for observed, apparent position, single time, and array time.
"""

from skyfield.api import load
from numpy import array
from tychos_skyfield import skyfieldlib as TS

eph_s = load('de421.bsp')

earth_ref = TS.ReferencePlanet('Earth', eph_s)
eph_t = TS.Ephemeris(earth_ref)
earth = eph_t["earth"]
jupiter = eph_t["jupiter"]

ts = load.timescale()

time = ts.tt(2020, 6, 21, 12, 0, 0)
radec_1 = earth.at(time).observe(jupiter).radec()
assert str(radec_1) == "(<Angle 19h 47m 50.19s>, <Angle -21deg 32' 13.2\">, <Distance 4.21064 au>)"
radec_2 = earth.at(time).observe(jupiter).apparent().radec()
assert str(radec_2) == "(<Angle 19h 47m 51.48s>, <Angle -21deg 32' 09.8\">, <Distance 4.21064 au>)"

times = ts.tt_jd(array([time.tt, time.tt + 100]))
radec_3 = earth.at(times).observe(jupiter).radec()
assert str(radec_3) == ("(<Angle 2 values from 19h 47m 50.19s to 19h 15m 43.97s>, "
                        "<Angle 2 values from -21deg 32' 13.2\" to -22deg 50' 25.3\">, "
                        "<Distance [4.21063884 4.81046205] au>)")
radec_4 = earth.at(times).observe(jupiter).apparent().radec()
assert str(radec_4) == ("(<Angle 2 values from 19h 47m 51.48s to 19h 15m 44.22s>, "
                        "<Angle 2 values from -21deg 32' 09.8\" to -22deg 50' 25.1\">, "
                        "<Distance [4.21063884 4.81046205] au>)")
