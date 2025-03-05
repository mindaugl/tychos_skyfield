"""
Test radec_direct calculation in Tychos coordinate system at default time 2000-06-21 12:00:00 and
 time of 1900-06-21 12:00:00, as calculated in 'date', 'june2000' and 'j2000June' epochs
"""

from tychos_skyfield import baselib as Tychos

system = Tychos.TychosSystem()
precision = 0.0001
string_format = "{}, {}, {:.4f}au"
mercury = system['Mercury']
earth = system['Earth']
polar_axis = system['polar_axis']

radec_date = mercury.radec_direct(earth, polar_axis, 'date')
assert string_format.format(*radec_date) == "7h 23m 42.48s, +20deg 53' 5.6\", 0.6285au"
radec_j2000June = mercury.radec_direct(earth, polar_axis, 'j2000June')
assert string_format.format(*radec_j2000June) == "7h 23m 42.48s, +20deg 53' 5.6\", 0.6285au"
radec_j2000 = mercury.radec_direct(earth, polar_axis, 'j2000')
assert string_format.format(*radec_j2000) == "7h 23m 40.78s, +20deg 53' 9.1\", 0.6285au"

system.move_system(2415192.0)
radec_date = mercury.radec_direct(earth, polar_axis, 'date')
assert string_format.format(*radec_date) == "7h 24m 15.32s, +22deg 29' 42.6\", 1.0903au"
radec_j2000June = mercury.radec_direct(earth, polar_axis, 'j2000June')
assert string_format.format(*radec_j2000June) == "7h 30m 19.73s, +22deg 16' 45.0\", 1.0903au"
radec_j2000 = mercury.radec_direct(earth, polar_axis, 'j2000')
assert string_format.format(*radec_j2000) == "7h 30m 18.01s, +22deg 16' 48.7\", 1.0903au"

system.move_system(2451545.0) # 2000-01-01 12:00:00 TT
assert(string_format.format(*mercury.radec_direct(earth, polar_axis, 'date')) ==
       string_format.format(*mercury.radec_direct(earth, None, 'j2000')))
