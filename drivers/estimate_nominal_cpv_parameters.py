from astropy import units as u, constants as c
import numpy as np

# Barnes09 gravity darkened model
# approximate shape of star as spheroid with oblateness gamma

# use representative params from Table1 of Gunther+22
# e.g. TIC 206544316, 45 myr

P = 7.73*u.hour
Teff = 3100*u.K
Rstar = 0.54*u.Rsun
Mstar = 0.22*u.Msun

P = 8*u.hour
Teff = 3000*u.K
Rstar = 0.5*u.Rsun
Mstar = 0.2*u.Msun

Ω = 2*np.pi / P

γ = Ω**2 * Rstar**3 / (2*c.G * Mstar)
R_cr = (c.G * Mstar / Ω**2)**(1/3)

R_l = c.c/Ω

g_eq = c.G * Mstar / (Rstar**2)
g_pole = c.G * Mstar / (((1+γ.cgs)*Rstar)**2)
beta = 0.25
T_pole = Teff * (g_pole/g_eq)**beta

v_eq = 2*np.pi*Rstar / P

print(f"Oblateness (a-b)/a ~= {γ.cgs}")
print(f"R_cr = {(R_cr/Rstar).cgs} Rstar")
print(f"R_cr = {(R_cr).to(u.AU)} AU")
print(f"R_cr = {(R_cr).to(u.Rsun)} Rsun")

print(f"R_l = {(R_l).to(u.Rsun)} Rsun")

print(f"g_eq/g_pole = {(g_eq/g_pole).cgs}")
print(T_pole.to(u.K))

print(f'v_eq = {v_eq.to(u.km/u.s)}')
