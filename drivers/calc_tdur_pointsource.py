import numpy as np
from astropy import units as u, constants as c

# a^3 = P^2 * GMtot / (4pi^2)

P = 18.56 * u.hr
R_star = 0.37 * u.Rsun
M_star = 0.22 * u.Msun
a = (P**2 * c.G * M_star / (4*np.pi**2))**(1/3)

Tdur = R_star * P / (np.pi * a)
print(f"{Tdur.to(u.hr):.2f}")

P = 18.56 * u.hr
R_star = (0.37+(0.1*0.37)) * u.Rsun
M_star = 0.22 * u.Msun
a = (P**2 * c.G * M_star / (4*np.pi**2))**(1/3)

Tdur = R_star * P / (np.pi * a)
print(f"{Tdur.to(u.hr):.2f}")
