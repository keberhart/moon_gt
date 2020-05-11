#!/usr/bin/env python
#
#   Trying to figure out if the lunar phase angle is incrreasing or decreasing
#   based on the moon phase description in skyfield
#
#------------------------------------------------------

from skyfield import almanac
from skyfield import api
from datetime import timedelta
from astropy import units as u

five_min = timedelta(days=2)

ts = api.load.timescale()
e = api.load('de421.bsp')

t0 = ts.now()

p0 = almanac.phase_angle(e, "moon", t0)
phase = almanac.moon_phases(e)

print p0.to(u.deg)
print almanac.MOON_PHASES[phase(t0)]
print "\n"

counter = 25
while counter > 0:

    t0 = ts.utc(t0.utc_datetime() + five_min)
    p0 = almanac.phase_angle(e, "moon", t0)
    
    print p0.to(u.deg)
    print almanac.MOON_PHASES[phase(t0)]
    print "\n"

    counter -= 1
