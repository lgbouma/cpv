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

P = 18.6*u.hour
Teff = 3000*u.K
Rstar = 0.5*u.Rsun
Mstar = 0.25*u.Msun

Ω = 2*np.pi / P

γ = Ω**2 * Rstar**3 / (2*c.G * Mstar)
R_cr = (c.G * Mstar / Ω**2)**(1/3)

R_l = c.c/Ω

g_eq = c.G * Mstar / (Rstar**2)
g_pole = c.G * Mstar / (((1+γ.cgs)*Rstar)**2)
logg = np.log10(g_eq.cgs.value)
beta = 0.25
T_pole = Teff * (g_pole/g_eq)**beta

v_eq = 2*np.pi*Rstar / P

T0 = P * Rstar / (R_cr * np.pi)
T0_stellarsurf = P * Rstar / (Rstar * np.pi) # at stellar_surface

Rp0 = 4*u.Rearth
δ0 = (Rp0/Rstar)**2
Rp1 = 2*u.Rearth
δ1 = (Rp1/Rstar)**2


print(f"P = {P}")
print(f"Oblateness (a-b)/a ~= {γ.cgs}")
print(f"R_cr = {(R_cr/Rstar).cgs} Rstar")
print(f"R_cr = {(R_cr).to(u.AU)}")
print(f"R_cr = {(R_cr).to(u.Rsun)} Rsun")

print(f"R_l = {(R_l).to(u.Rsun)} Rsun")

print(f"g_eq/g_pole = {(g_eq/g_pole).cgs}")
print(T_pole.to(u.K))

print(f"logg = {logg:.4f}")

print(f'v_eq = {v_eq.to(u.km/u.s)}')

print(f'T0 = {T0.to(u.hr)}')
print(f'T0 [phase units] = {(T0.to(u.hr)/P.to(u.hr)).cgs.value:.4f}')

print(f'Assuming Rstar = {Rstar}')
print(f"For Rp0={Rp0}, δ={δ0.cgs:.2e}")
print(f"For Rp1={Rp1}, δ={δ1.cgs:.2e}")
