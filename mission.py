import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris

from poliastro.bodies import Sun, Earth, Jupiter
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver
from poliastro.iod import izzo
from poliastro.plotting import plot, OrbitPlotter
from poliastro.util import norm

solar_system_ephemeris.set("jpl")

## Initial data
# Links and sources: https://github.com/poliastro/poliastro/wiki/EuroPython:-Per-Python-ad-Astra
date_launch = Time("2011-08-05 16:25", scale='utc')
C_3 = 31.1 * u.km**2 / u.s**2 #specific orbital energy when orbiting Earth's surface (en.wikipedia.org/wiki/Specific_Orbital_Energy)
date_flyby = Time("2013-10-09 19:21", scale='utc')
date_arrival = Time("2016-07-05 03:18", scale='utc')

# Initial state of the Earth
ss_e0 = Orbit.from_body_ephem(Earth, date_launch)
r_e0, v_e0 = ss_e0.rv()

#r_e0
#v_e0

# State of the Earth the day of the flyby
ss_efly = Orbit.from_body_ephem(Earth, date_flyby)
r_efly, v_efly = ss_efly.rv()

# Assume that the insertion velocity is tangential to that of the Earth
dv = C_3**.5 * v_e0 / norm(v_e0) #project sqrt(C_3) onto unit vector tangential to Earth for change in velocity
man = Maneuver.impulse(dv) #set Poliastro "maneuver" to the change in velocity dv

# Inner Cruise 1
ic1 = ss_e0.apply_maneuver(man) #set inner orbit to leave initial state in Earth's orbit with defined "maneuver" at each step
ic1.rv()

ic1.period.to(u.year) #convert period of inner cruise 1 orbit in years

op = OrbitPlotter()

op.plot(ss_e0) #plot Earth's orbit
op.plot(ic1) #plot inner cruise 1 orbit

# We propagate until the aphelion
ss_aph = ic1.propagate(ic1.period / 2) #aphelion: the point in the orbit at which it is furthest from the sun = ic1.period/2
ss_aph.epoch #return time of aphelion

# Let's compute the Lambert solution to do the flyby of the Earth
# Lambert's problem is concerned with the determination of an orbit from two position vectors and the time of flight
time_of_flight = date_flyby - ss_aph.epoch # time of aphelion to flyby of Earth
time_of_flight

# Determine the orbit from gravitational constant of sun, position of aphelion, position of flyby, and time of flight
(v_aph, v_fly), = izzo.lambert(Sun.k, ss_aph.r, ss_efly.r, time_of_flight) 

# Check the difference between initial velocity necessary for orbit and velocity at aphelion
norm(v_aph - ss_aph.v)  # Difference is too great...intermediate orbit defined below

ss_aph_post = Orbit.from_vectors(Sun, ss_aph.r, v_aph, epoch = ss_aph.epoch)
ss_junofly = Orbit.from_vectors(Sun, r_efly, v_fly, epoch = date_flyby)

op = OrbitPlotter()

op.plot(ss_e0, label="Earth")
op.plot(ic1, label="Inner Cruise 1")
#op.plot(ss_efly)
op.plot(ss_aph_post, label="Back to Earth")

# And now, go to Jupiter!
ss_j = Orbit.from_body_ephem(Jupiter, date_arrival)
r_j, v_j = ss_j.rv()

(v_flypre, v_oip), = izzo.lambert(Sun.k, r_efly, r_j, date_arrival - date_flyby)

ss_oip = Orbit.from_vectors(Sun, r_j, v_oip, epoch=date_flyby)

fig, ax = plt.subplots(figsize=(9, 12))

op = OrbitPlotter(ax)

op.plot(ss_e0, label="Earth")
op.plot(ic1, label="Inner Cruise 1")
#op.plot(ss_efly)
op.plot(ss_aph_post, label="Back to Earth")
op.plot(ss_oip, label="Jupiter Orbit Insertion Phase")
op.plot(ss_j, label="Jupiter")

fig.savefig("jupiter.png")
plt.show()
