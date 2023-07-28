from astropy import units as u, constants as c
import numpy as np

# way #1: set magnetic pressure B^2/8pi equal to 3/2 nkT

B_0 = 1e3*u.cm**(-1/2)*u.g**(1/2)*u.s**(-1)

T = 1e6*u.K
n = 1e12/(u.cm**3)  # !!!! this is high!?   [seems super high!]
m = n*c.m_p
print(m.cgs)

R_star = 0.4*u.Rsun

r_A = R_star * (
    B_0**2 / (3*np.pi) * (1/(n*c.k_B*T))
).cgs**(1/4)

print(f'r_A = {r_A:.4f}')
print(f'r_A/r_* = {r_A/R_star:.4f}')

# way #2: hartmann+2016
B_3 = 1e4 / (1e3)
R_2 = 0.37 / 2
M_0pt5 = 0.22 / 0.5
Mdot_m8 = 1e-8 / 1e-8 # outflow rate.  ryden ch11 "the solar wind" asserts 1e-7 - 1e-9 for TT.  maybe 1e-10 or 1e-11 allowed for wTTs

r_M = 18 * (
    B_3 ** (4/7) *
    R_2 ** (12/7) *
    M_0pt5 ** (1/7) *
    Mdot_m8 ** (2/7)
)*u.Rsun

print(f'r_M = {r_M:.2f}')
print(f'r_M/r_* = {r_M/R_star:.2f}')

R_cr_by_Rstar = 5.76 # tic 4029*
print(f'r_CR/r_* = {R_cr_by_Rstar}')

print(f'r_M > r_CR (centrifugal magnetosphere) = {r_M/R_star > R_cr_by_Rstar}')

# same but for a rapidly rotating 1 Msun star
# RESULT: you either need to tune down B, or lower Mdot.  Because the radius
# scaling goes the _wrong_ way to help you!   Kochukhov2021 suggests that
# tuning down B is plausible.
# ... could also just be a lifetime and IMF thing.  If you look at a large
# enough sample of PMS G dwarfs, they should meet this condition too, I think.
# Radii are big!!  And that's the largest dependence.
B_3 = (300) / (1e3)  # Kochukhov "typical"
R_2 = 1 / 2
M_0pt5 = 1 / 0.5
Mdot_m8 = 1e-9 / 1e-8 # outflow rate.  ryden ch11 "the solar wind" asserts 1e-7 - 1e-9 for TT.  maybe 1e-10 or 1e-11 allowed for wTTs

r_M = 18 * (
    B_3 ** (4/7) *
    R_2 ** (12/7) *
    M_0pt5 ** (1/7) *
    Mdot_m8 ** (2/7)
)*u.Rsun

print(f'r_M = {r_M:.2f}')
print(f'r_M/r_* = {r_M/R_star:.2f}')

R_cr_by_Rstar = 5.76 # tic 4029*
print(f'r_CR/r_* = {R_cr_by_Rstar}')

print(f'r_M > r_CR (centrifugal magnetosphere) = {r_M/R_star > R_cr_by_Rstar}')

