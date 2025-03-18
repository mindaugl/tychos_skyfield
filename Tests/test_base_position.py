"""
Test location in Tychos coordinate system for the default time of 2000-06-21 12:00:00
"""

from tychos_skyfield import baselib as Tychos

system = Tychos.TychosSystem()
precision = 0.0001

location = system['Jupiter'].location
assert abs(location[0] - 296.65674) < precision
assert abs(location[1] - (-11.07307)) < precision
assert abs(location[2] - (-515.52472)) < precision

location = system['Mars'].location
assert abs(location[0] - (-16.02508)) < precision
assert abs(location[1] - 2.68133) < precision
assert abs(location[2] - (-297.02498)) < precision
