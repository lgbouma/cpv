import numpy as np
import astropy.units as u
import astropy.constants as const

def calculate_cpv_mass(star_radius, dip_depth, particle_radius, particle_density):
    # Convert inputs to cgs units
    R_star = star_radius.to(u.cm)
    r_particle = particle_radius.to(u.cm)
    rho_particle = particle_density.to(u.g / u.cm**3)

    # Calculate area of obscuring material
    area_obscuring = dip_depth * np.pi * R_star**2

    # Calculate optical properties
    sigma = np.pi * r_particle**2
    L = np.sqrt(area_obscuring)

    # Calculate number density for optical depth >= 1
    n = 1 / (sigma * L)
    print(f"Number density n: {n:.2e}")

    # Calculate total number of particles
    N = n * area_obscuring * L

    # Calculate mass of each particle
    m_particle = 4/3 * np.pi * r_particle**3 * rho_particle

    # Calculate total mass
    M_total = N * m_particle

    return M_total

# Set up the parameters
R_star = 0.3 * const.R_sun
dip_depth = 0.01
r_particle = 1 * u.micron
rho_particle = 3 * u.g / u.cm**3  # Now in cgs units

# Calculate the mass
M = calculate_cpv_mass(R_star, dip_depth, r_particle, rho_particle)

# Print results
print(f"Total mass of obscuring material: {M:.2e}")
print(f"In Earth masses: {(M / const.M_earth).decompose():.2e}")
print(f"Fraction of Ceres mass: {(M / (9.3e20 * u.kg)).decompose():.2e}")
