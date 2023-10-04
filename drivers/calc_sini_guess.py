import numpy as np
from astropy import constants as c, units as u

# tic 89026133 estimate
Rstar = 1.33*u.Rsun
P = 11.2*u.hr
v = (2*np.pi*Rstar / P).to(u.km/u.s)
vsini = 60 *u.km/u.s

print(v)
print(vsini)
print(vsini/v)

# tic 167664935 estimate
Rstar = 1.42*u.Rsun
P = 14.05*u.hr
v = (2*np.pi*Rstar / P).to(u.km/u.s)
vsini = 100 *u.km/u.s

# tic 141146667
Rstar = 0.42*u.Rsun
P = 3.93*u.hr
v = (2*np.pi*Rstar / P).to(u.km/u.s)

print(v)
print(vsini)
print(vsini/v)
