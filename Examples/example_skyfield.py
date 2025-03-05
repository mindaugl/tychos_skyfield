"""
An example script showing how to use Tychos python model with Skyfield, where tychos objects
can be treated as planets in Skyfield (JPL ephemeris data).
Documentation of skyfield functionality can be found at: https://rhodesmill.org/skyfield/
NOTE: Ra/Dec calculations in this example are done in the ICRF/J2000 frame. Skyfield supports
other frames too, including 'date' frame co-moving with earth, that is used in the Tychosium,
but the Skyfield uses Earth 'nutation' model to transform to the 'date' frame which introduces
about 15 arcseconds discrepancy as compared to the Tychosium calculation.
"""

from skyfield.api import load
from tychos_skyfield import skyfieldlib as TS

# Get Tychos observable objects:
print("Tychos observable objects:", TS.TychosSkyfield.get_observable_objects())

# Load skyfield planets data:
skyfield_objs = load('de421.bsp')

# Create Earth Skyfield object and reference object for Tychos:
earth_s = skyfield_objs['Earth']
earth_ref = TS.ReferencePlanet('Earth', skyfield_objs)

# Tychos Jupiter object (with Earth as reference object) that complies with Skyfield infrastructure:
jupiter_t = TS.TychosSkyfield("Jupiter", earth_ref)

# Set Terrestrial Time (corresponds to the UTC time used in the Tychos model, same Julian Day):
ts = load.timescale()
time = ts.tt(2020, 6, 21, 12, 0, 0)

# Observe position and Ra/Dec/Dist in the ICRF/J2000 frame as observed from Earth:
print("Instantaneous Position:", jupiter_t.at(time).position)
print("Instantaneous Ra/Dec/Dist:       ", jupiter_t.at(time).radec())
print("Observed Ra/Dec/Dist:            ", earth_s.at(time).observe(jupiter_t).radec())
print("Observed Apparent Ra/Dec/Dist:   ", earth_s.at(time).observe(jupiter_t).apparent().radec())

# Check that instantaneous Ra/Dec/Dist agree with native tychos calculation in the ICRF/J2000 frame:
jupiter_t.move_system(time.tt)
print("\nInstantaneous Ra/Dec/Dist native:",
      jupiter_t.native_object().radec_direct(jupiter_t.native_object("Earth")))

# Observe Skyfield (JPL ephemeris) Jupiter in the ICFR/J2000 frame:
print("\nObserved RA/Dec/Dist skyfield:   ",
      earth_s.at(time).observe(skyfield_objs["Jupiter_barycenter"]).radec())
