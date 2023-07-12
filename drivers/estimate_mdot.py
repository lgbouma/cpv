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

# What about for growth seen in TIC 224283342, of 3% in depth over 20 days?
# from Sanderson2023 eq23
τ = 1
Qext = 3
rdust = 1 # micron
rhodust = 3000 # kg m^-3

delta = 0.03
rcloud_by_rstar = np.sqrt(delta)

Rstar = 0.39*u.Rsun

m_cloud = (
    1.6e12 * (τ/1) * (Qext/3)**(-1)
    * ( (1/0.1) * (rcloud_by_rstar) * (Rstar/ (0.34*u.Rsun) ) )**2
    * (rdust/1) * (rhodust / 3000)
)*u.kg
print(f"Mcloud {m_cloud:.2e}")

dt = 20*u.day

dM_dt = m_cloud / dt
print(f'Mdot {dM_dt.to(u.Mearth/u.yr):.1e}')
