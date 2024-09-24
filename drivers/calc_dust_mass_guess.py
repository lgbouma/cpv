import numpy as np
from astropy import units as u
from astropy.constants import sigma_sb

# Constants
T_star = 3000 * u.K  # Temperature of the star in K
T_dust = 1500 * u.K  # Temperature of the dust in K
R_star = 0.4 * u.Rsun  # Radius of the star in solar radii
a = 1 * u.micron  # Radius of the dust grains in microns
rho = 3 * (u.g / u.cm**3)  # Density of silicate grains in g/cm^3
Q_abs = 1  # Absorption efficiency (assumed to be 1 for simplicity)
Q_em = 1  # Emission efficiency (assumed to be 1 for simplicity)

# Calculate the area of the annulus
R_in = 2.9 * R_star
R_out = 3.2 * R_star
A_annulus = np.pi * (R_out**2 - R_in**2)

# Calculate the luminosity of the star
L_star = 4 * np.pi * R_star**2 * sigma_sb * T_star**4

# Calculate the luminosity of the dust
L_dust = A_annulus * sigma_sb * T_dust**4

# Calculate the ratio of luminosities
L_ratio = (T_dust / T_star)**4 * (A_annulus / (4 * np.pi * R_star**2))

# Calculate the cross-sectional area of a single dust grain
A_dust_grain = np.pi * a**2

# Calculate the mass of a single dust grain
V_dust_grain = (4/3) * np.pi * a**3
m_dust_grain = V_dust_grain * rho

# Calculate the number of dust grains needed
N_dust_grains = L_dust / (A_dust_grain * sigma_sb * T_dust**4)

# Calculate the total mass of the dust
M_dust = N_dust_grains * m_dust_grain

IN_ABSORPTION = 0
if IN_ABSORPTION:
    # Calculate the mass absorption coefficient
    # kappa: typically cm^2 / g - here's it's 2500 cm^2/g
    kappa = (3 * Q_abs) / (4 * a * rho)

    # Calculate the mass of the dust using nutty dimensional scaling argument.
    M_dust = A_annulus / kappa

# Display results with correct units
print(
    f"R_star = {R_star:.2f}\n"
    f"R_in = {R_in:.2f}\n"
    f"R_out = {R_out:.2f}\n"
    f"T_dust = {T_dust}\n"
    f"L_star = {L_star.to(u.Lsun):.1e},\n"
    f"L_dust = {L_dust.to(u.Lsun):.1e},\n"
    f"L_ratio = {L_ratio:.3f},\n"
    f"N_dust_grains = {N_dust_grains.cgs:.1e}\n"
    f"M_dust = {M_dust.to(u.kg):.1e},\n"
    f"M_dust = {M_dust.to(u.Mearth):.1e}\n"
)

