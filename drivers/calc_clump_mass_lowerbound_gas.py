"""
Assuming κ=0.1cm^2/g (thompson scattering opacity == electron scattering
opacity, assuming hot fully ionized plasma).

Actual κ_thompson = σ_Thomspon / m_proton = 0.4 cm^2 / g

However this is reduced by 2-4x in practical applications b/c of heavier
elements and partial ionization.
"""
import numpy as np
import astropy.units as u
import astropy.constants as const

def calculate_cpv_gas_mass(star_radius, dip_depth, gas_opacity):
    # Convert inputs to cgs units
    R_star = star_radius.to(u.cm)
    
    # Calculate area of obscuring material
    area_obscuring = dip_depth * np.pi * R_star**2
    
    # Estimate path length
    L = np.sqrt(area_obscuring)
    
    # Calculate required density for optical depth >= 1
    rho = 1 / (gas_opacity * L)
    
    # Calculate total mass
    M_total = rho * area_obscuring * L
    
    return M_total

# Set up the parameters
R_star = 0.3 * const.R_sun
dip_depth = 0.01
gas_opacity = 0.1 * u.cm**2 / u.g  # Opacity for hot, ionized gas

# Calculate the mass
M = calculate_cpv_gas_mass(R_star, dip_depth, gas_opacity)

# Print results
print(f"Total mass of obscuring gas: {M:.2e}")
print(f"In Earth masses: {(M / const.M_earth).decompose():.2e}")
print(f"Fraction of Ceres mass: {(M / (9.3e20 * u.kg)).decompose():.2e}")
