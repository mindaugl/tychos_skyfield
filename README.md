# Welcome to the tychos_skyfield repository

Tychos (Tychosium) model implementation in Python3 and integration with the Skyfield package.

This package enables users to run Tychos solar system model in Python and to import it into Skyfield package.

## Tychos model

```
  from tychos_skyfield import baselib as T
    
  # Create solar system:
  system = T.TychosSystem()
  
  # Move system to particular time (2020-06-21 12:00:00):
  julian_day = 2459022.0
  system.move_system(julian_day)
  
  # Check list of observable objects:
  print("Observable objects:", system.get_observable_objects())
  
  # Check location of object, given in centi-AU (earth-sun distance is ~100):
  print("Jupiter location:", system['Jupiter'].location)
  
  # Get RA/Dec/Distance as calculated in Tychosium (in the frame associated with the Earth's orientation):
  ra, dec, dist = system['Jupiter'].radec_direct(system['Earth'], system['Polar_axis'], 'date')
  print("RA/Dec/Distance in the moving frame:", ra, ",", dec, ",", dist)
  
  # Get RA/Dec/Distance for the J2000 epoch (in the frame of Earth being on the date 2000-01-01):
  ra, dec, dist = system['Jupiter'].radec_direct(system['Earth'], None, 'j2000')
  print("RA/Dec/Distance in the j2000 frame: ", ra, ",", dec, ",", dist)
```

With result:
```
  Observable objects: ['sun', 'mercury', 'moon', 'venus', 'mars', 'phobos', 'deimos', 'jupiter', 'saturn', 'uranus', 'neptune', 'halleys', 'eros']
  Jupiter location: [177.79634343  -3.76073738 343.90806095]
  RA/Dec/Distance in the moving frame: 19h 49m 4.52s , -21deg 28' 57.8" , 4.210632538170995
  RA/Dec/Distance in the j2000 frame:  19h 47m 50.81s , -21deg 32' 11.8" , 4.210632538170995
```

## Tychos integration with Skyfield

```
  from skyfield.api import load
  from tychos_skyfield import skyfieldlib as TS
  
  # Get Tychos observable objects:
  print("Tychos observable objects:", TS.TychosSkyfield.get_observable_objects())
  
  # Load skyfield planets data (~17 MB):
  skyfield_objs = load('de421.bsp')
  
  # Create Earth Skyfield object and reference object for Tychos:
  earth_s = skyfield_objs['Earth']
  earth_ref = TS.ReferencePlanet('Earth', skyfield_objs)
  
  # Tychos Jupiter object (with Earth as reference object) that complies with Skyfield infrastructure:
  jupiter_t = TS.TychosSkyfield("Jupiter", earth_ref)
  
  # Set time:
  ts = load.timescale()
  time = ts.tt(2020, 6, 21, 12, 0, 0)
  
  # Observe position and Ra/Dec/Dist in the ICRF/J2000 frame as observed from Earth:
  print("Jupiter instantaneous Position:", jupiter_t.at(time).position)
  print("Instantaneous Ra/Dec/Dist:       ", jupiter_t.at(time).radec())
  print("Observed Ra/Dec/Dist:            ", earth_s.at(time).observe(jupiter_t).radec())
  print("Observed Apparent Ra/Dec/Dist:   ", earth_s.at(time).observe(jupiter_t).apparent().radec())
  
  # Check that instantaneous Ra/Dec/Dist agree with native tychos calculation in the ICRF/J2000 frame:
  jupiter_t.move_system(time.tt)
  print("Instantaneous Ra/Dec/Dist native:", jupiter_t.native_object().radec_direct(jupiter_t.native_object("Earth")))
```
With result:
```
  Jupiter instantaneous Position: [ 1.77579358 -3.49095666 -1.54570538] au
  Instantaneous Ra/Dec/Dist:        (<Angle 19h 47m 50.81s>, <Angle -21deg 32' 11.8">, <Distance 4.21063 au>)
  Observed Ra/Dec/Dist:             (<Angle 19h 47m 50.19s>, <Angle -21deg 32' 13.2">, <Distance 4.21064 au>)
  Observed Apparent Ra/Dec/Dist:    (<Angle 19h 47m 51.48s>, <Angle -21deg 32' 09.8">, <Distance 4.21064 au>)
  Instantaneous Ra/Dec/Dist native: ('19h 47m 50.81s', '-21deg 32\' 11.8"', np.float64(4.210632538170995))
```

## Installation

```
pip install tychos_skyfield
```

The dependencies are:
-  numpy
-  scipy
-  skyfield (only for tychos skyfield library)

The package was tested to work with python 3.12.3, numpy 2.2.1, scipy 1.15.1, skyfield 1.51.

## References

Tychosium 3-D (source code): https://codepen.io/pholmq/full/XGPrPd

Tychos project: https://www.tychos.space/

Skyfield: https://rhodesmill.org/skyfield/
