from astropy import constants as c, units as u
import numpy as np

# Mstar, Rstar, Pstar,
params = [
    [0.3*u.Msun,  0.4*u.Rsun, 1*u.day],
    [0.3*u.Msun,  0.4*u.Rsun, 0.2*u.day],
    [0.22*u.Msun,  0.54*u.Rsun, 7.7*u.hour], # TIC 206544316
    [0.10*u.Msun, 0.27*u.Rsun, 3.6*u.hour], # TIC 201789285
]

for p in params:

    M_star = p[0]
    R_star = p[1]
    P_star = p[2]

    a_corotation = (P_star/R_star)**2 * (c.G * M_star / (4*np.pi**2))

    msg = (
        f"M_star={M_star}, R_star={R_star}, P_star={P_star.to(u.day):.2f}, "+
        f"a_corotation={a_corotation.to(u.Rsun):.2f}"
    )

    print(msg)
