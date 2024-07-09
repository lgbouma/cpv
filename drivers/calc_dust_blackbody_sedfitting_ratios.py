#!/usr/bin/env python3
"""
calc_dust_mass_upperbound_sedfitting.py

This script calculates the maximum luminosity of a Teff=1000 K blackbody that could be 
fit under the curve of a Teff=3000 K blackbody while meeting observational constraints.
It generates a log-log plot of flux vs. wavelength to show the result.

Usage:
    python calc_dust_mass_upperbound_sedfitting.py
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling.physical_models import BlackBody
from astropy import units as u, constants as c

def blackbody_flux(temperature, wavelengths):
    """
    Calculate the blackbody flux for a given temperature and array of wavelengths.

    Args:
        temperature (float): Temperature in Kelvin.
        wavelengths (numpy.ndarray): Wavelengths in meters.

    Returns:
        numpy.ndarray: Blackbody flux B_nu(T) in (erg / (cm ** 2 * s * Hz * sr), with
        scale=1, and also F_bol (in erg / cm^2 s).
    """
    bb = BlackBody(temperature=temperature * u.K)

    B_nu = bb(wavelengths)

    bolometric_flux = bb.bolometric_flux

    return B_nu, bolometric_flux


def calculate_flux_ratio(teff_star, teff_dust):
    """
    Calculate the flux ratio of two blackbodies.

    Args:
        teff_star (float): Effective temperature of the star in Kelvin.
        teff_dust (float): Effective temperature of the dust in Kelvin.

    Returns:
        float: Maximum luminosity ratio of the dust blackbody to the star blackbody.
    """
    wavelengths = np.logspace(2, 5, 1000) * u.angstrom

    flux_star, flux_star_bol = blackbody_flux(teff_star, wavelengths)
    flux_dust, flux_dust_bol = blackbody_flux(teff_dust, wavelengths)

    # scale it at ~3 microns, where the dust brightness peaks.  (since we
    # assumed Teff=1500K).
    ind = np.argmax(flux_dust)

    ratio = flux_star[ind] / flux_dust[ind]

    dust_bb_scaling_factor = 0.05*ratio.cgs.value
    print(dust_bb_scaling_factor)

    return dust_bb_scaling_factor

def plot_blackbodies(teff_star, teff_dust, dust_bb_scaling_factor, output_file):
    """
    Plot the blackbody curves for the star and the dust.

    Args:
        teff_star (float): Effective temperature of the star in Kelvin.
        teff_dust (float): Effective temperature of the dust in Kelvin.
        dust_bb_scaling_factor (float): Scaling factor for the dust blackbody luminosity.
        output_file (str): Path to save the plot.
    """
    wavelengths = np.logspace(2, 5, 1000) * u.angstrom

    flux_star, flux_star_bol = blackbody_flux(teff_star, wavelengths)
    flux_dust, flux_dust_bol = blackbody_flux(teff_dust, wavelengths)
    flux_dust *= dust_bb_scaling_factor
    flux_dust_bol *= dust_bb_scaling_factor

    plt.figure(figsize=(10, 6))
    plt.loglog(wavelengths, flux_star, label=f'Star (Teff={teff_star} K)')
    plt.loglog(wavelengths, flux_dust, label=f'Dust (Teff={teff_dust} K, scaled {dust_bb_scaling_factor:.1e})')
    flux_total = flux_star + flux_dust
    plt.loglog(wavelengths, flux_total, label=f'Total', c='k')
    plt.xlabel('Wavelength (A)')
    plt.ylabel('F_nu (odd units)')
    ymax = flux_total.max()*1.2
    ymin = 0.01*ymax
    plt.ylim((ymin.value, ymax.value))
    plt.xlim((4000, 1e5))
    plt.legend()
    plt.title('Blackbody Spectra')
    plt.grid(True, which="both", ls="--")
    plt.savefig(output_file)

def calc_flux_ratio_at_2microns(teff_star, teff_dust, dust_bb_scaling_factor):
    """
    Calculate the flux ratio at 2 microns (20000 angstroms).

    Args:
        teff_star (float): Effective temperature of the star in Kelvin.
        teff_dust (float): Effective temperature of the dust in Kelvin.
        dust_bb_scaling_factor (float): Scaling factor for the dust blackbody luminosity.

    Returns:
        float: Flux ratio at 2 microns.
    """
    wavelength_2microns = 2e-6 * u.m
    flux_star, _ = blackbody_flux(teff_star, wavelength_2microns)
    flux_dust, _ = blackbody_flux(teff_dust, wavelength_2microns)
    flux_dust *= dust_bb_scaling_factor
    flux_total = flux_star + flux_dust
    return (flux_star / flux_total).value

def bol_flux_ratio(teff_star, teff_dust, dust_bb_scaling_factor):
    """
    Calculate the bolometric flux ratio of the scaled 1500 K blackbody to the 3000 K blackbody.

    Args:
        teff_star (float): Effective temperature of the star in Kelvin.
        teff_dust (float): Effective temperature of the dust in Kelvin.
        dust_bb_scaling_factor (float): Scaling factor for the dust blackbody luminosity.

    Returns:
        float: Luminosity ratio of the scaled 1500 K blackbody to the 3000 K blackbody.
    """
    wavelengths = np.logspace(2, 5, 1000) * u.angstrom

    flux_star, flux_star_bol = blackbody_flux(teff_star, wavelengths)
    flux_dust, flux_dust_bol = blackbody_flux(teff_dust, wavelengths)
    flux_dust *= dust_bb_scaling_factor
    flux_dust_bol *= dust_bb_scaling_factor

    dust_to_star_bol_flux_ratio = flux_dust_bol / flux_star_bol

    return dust_to_star_bol_flux_ratio

def calculate_dust_mass(L_dust, T_dust, dust_bb_scaling_factor, wavelength=2e-6):
    """
    Calculate the dust mass given the dust luminosity and temperature.

    Args:
        L_dust (float): Dust luminosity in solar luminosities.
        T_dust (float): Dust temperature in Kelvin.
        wavelength (float): Wavelength in meters for κ_ν calculation. Default
        is 2 micron.

    Returns:
        float: Dust mass in kg.
    """

    # Calculate κ_ν for 1 micron silicate grains (approximate value)
    # This is a rough approximation and may need adjustment based on specific dust properties
    # 
    # The typical values for κ_ν's prefactor (κ_0) in cgs units can vary depending on
    # the specific dust composition, grain size distribution, and wavelength regime.
    # 
    # For interstellar dust in the submillimeter to millimeter wavelength range:
    # κ_0 is typically in the range of 0.1 - 10 cm²/g at λ_0 = 850 μm or 1 mm.
    # For protoplanetary disks:
    # A commonly used value is κ_0 ≈ 2 - 3 cm²/g at λ_0 = 850 μm or 1 mm.
    # For debris disks:
    # Values can be higher, with κ_0 ≈ 1 - 10 cm²/g at λ_0 = 850 μm.
    # For silicate-dominated dust:
    # At shorter wavelengths (e.g., λ_0 = 10 μm), κ_0 can be around 1000 cm²/g.
    # For carbonaceous grains:
    # These can have higher opacities, with κ_0 potentially reaching 100 - 1000 cm²/g
    # in the far-infrared.
    # 
    # In many astrophysical studies, researchers often use:
    # κ_ν = κ_0 (ν / ν_0)^β
    # where β is the dust opacity index, typically between 1 and 2 for interstellar dust.
    # For your calculations involving 1 micron grains, you might want to use a κ_0
    # value on the higher end of these ranges, as smaller grains generally have
    # higher opacities. However, the exact value would depend on the composition and
    # the specific wavelength you're most interested in.

    kappa_nu = 1e3 * (u.cm**2/u.g) * (wavelength / 1e-6)**(-1.5)

    # Calculate B_ν(T_dust) using the blackbody_flux function
    B_nu = dust_bb_scaling_factor * blackbody_flux(T_dust, np.array([wavelength]))[0]

    # Calculate dust mass
    M_dust = (
        L_dust.cgs.value / (4 * np.pi * kappa_nu.cgs.value * B_nu.cgs.value)
    ) * u.g

    print(M_dust.to(u.Msun))
    import IPython; IPython.embed()

    return M_dust.to(u.Msun)

def main():
    teff_star = 3000
    teff_dust = 1500
    L_star = 0.01 * u.L_sun
    output_file = '../results/dust_oom_estimates/blackbody_spectra.png'

    dust_bb_scaling_factor = calculate_flux_ratio(teff_star, teff_dust)
    plot_blackbodies(teff_star, teff_dust, dust_bb_scaling_factor, output_file)

    flux_ratio_at_2microns = calc_flux_ratio_at_2microns(
        teff_star, teff_dust, dust_bb_scaling_factor
    )
    print(f"Flux ratio at 2 microns (3000 K / total): {flux_ratio_at_2microns:.4f}")

    dust_to_star_bol_flux_ratio = bol_flux_ratio(teff_star, teff_dust, dust_bb_scaling_factor)
    print(f"Bol flux ratio (scaled 1500 K / 3000 K): "
          f"{dust_to_star_bol_flux_ratio:.4f}")


if __name__ == '__main__':
    main()
