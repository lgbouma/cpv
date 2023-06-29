from astropy import units as u, constants as c
import numpy as np

# imagine a state-switch represents a mass loss event.
# ~3 seen over 6 months in LP 12-502.
dN_dt = 3 / (0.5*u.year)

# Sanderson+2023
M = 1e12 * u.kg

# Pitjeva2018, found via wikipedia(!)
# 3% of the moon's mass.
# For comparison, the AU Mic debris disk, spread over 200AU, is thoguht to
# contain at least a lunar mass worth of material.
# https://ui.adsabs.harvard.edu/abs/2005ApJ...634.1372C
M_asteroidbelt = 2.4e21 * u.kg

Mdot = M * dN_dt

print(f'Mdot {Mdot.to(u.Msun/u.yr):.1e}')
print(f'Mdot {Mdot.to(u.Mearth/u.yr):.1e}')

timescale = 1e8 * u.yr

Mcum = Mdot * timescale

print(f'Mcum over 1e8 yr: {Mcum.to(u.Mearth):.1e}')
print(f'Mcum over 1e8 yr divided by asteroid belt mass: {Mcum/M_asteroidbelt:.2f}')
