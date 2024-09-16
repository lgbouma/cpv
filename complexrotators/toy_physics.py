"""
Contents:

Assuming dust cloud:
    dust_transit_depth

Assuming plasma cloud + transit light source effect:
    planck_lambda
    stellar_flux
    transit_depth_with_cools
"""
import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy.constants import h, c, k

def dust_transit_depth(
    wavelength: np.ndarray,
    beta: float,
    r_star: float,
    Sigma: float,
    A_clump_ratio: float
    ) -> np.ndarray:
    """
    Calculates the transit depth due to an occulting dust cloud as a function
    of wavelength assuming that

    δ(λ) = (A_clump / A_star) * (1 - exp[- κ(λ) * Σ])

    Args:
        wavelength (np.ndarray): Array of wavelengths (in microns).
        beta (float): Power-law index for the opacity (1 <= beta <= 2).
        r_star (float): Radius of the star (in solar radii).
        Sigma (float): Surface density of the clump.
        A_clump_ratio (float): Ratio of clump area to star area.

    Returns:
        np.ndarray: Dip depth at each wavelength.
    """
    # Star area
    A_star = np.pi * (r_star * 6.96e8)**2  # Convert solar radii to meters

    # Opacity power-law, assume kappa(λ) ∝ λ^(-beta)
    kappa_1micron = 1e3  # cm^2/g
    kappa = kappa_1micron * (wavelength / 1.0)**(-beta)  # Physical normalization

    # Optical depth τ(λ) = kappa(λ) * Sigma
    tau = kappa * Sigma

    # Blocked area is proportional to the clump area and depends on optical depth
    A_clump = A_clump_ratio * A_star
    A_blocked = A_clump * (1 - np.exp(-tau))

    # Dip depth δ(λ) = A_blocked / A_star
    delta = A_blocked / A_star

    return delta

def planck_lambda(T: float, wavelength: np.ndarray) -> np.ndarray:
    """
    Planck's law to calculate the blackbody spectral radiance at a given
    temperature and wavelength.

    Args:
        T (float): Temperature in Kelvin.
        wavelength (np.ndarray): Wavelength in microns.

    Returns:
        np.ndarray: Blackbody flux at each wavelength.
    """
    # Convert wavelength to meters
    wavelength_m = wavelength * 1e-6

    # Planck's law formula
    B_lambda = (2.0 * h * c**2) / (wavelength_m**5) / (np.exp((h * c) / (wavelength_m * k * T)) - 1)

    return B_lambda

def stellar_flux(f: float, T_hot: float, T_cool: float, wavelength: np.ndarray) -> np.ndarray:
    """
    Calculate the total stellar flux accounting for cool regions and hot regions.

    Args:
        f (float): Fraction of the stellar surface covered by cool regions.
        T_hot (float): Temperature of the hot regions in Kelvin.
        T_cool (float): Temperature of the cool regions in Kelvin.
        wavelength (np.ndarray): Wavelength in microns.

    Returns:
        np.ndarray: Total stellar flux at each wavelength.
    """
    B_hot = planck_lambda(T_hot, wavelength)
    B_cool = planck_lambda(T_cool, wavelength)

    # Total flux with fractional contributions from cool regions and hot regions
    F_star = (1 - f) * B_hot + f * B_cool

    return F_star

def transit_depth_tlse(
    f: float,
    T_hot: float,
    T_cool: float,
    delta_gray: float,
    wavelength: np.ndarray
    ) -> np.ndarray:
    """
    Calculate the transit depth as a function of wavelength, accounting for the
    transit light source effect (cool regions + hot regions, with a gray
    occulter).

    Args:
        f (float): Fraction of the stellar surface covered by cool regions.
        T_hot (float): Temperature of the hot regions in Kelvin.
        T_cool (float): Temperature of the cool regions in Kelvin.
        delta_gray (float): The gray transit depth.
        wavelength (np.ndarray): Wavelength in microns.

    Returns:
        np.ndarray: Dip depth at each wavelength.
    """
    # Stellar flux as seen from the entire star (hot regions + cool regions)
    F_star = stellar_flux(f, T_hot, T_cool, wavelength)

    # Flux from the transit chord (hot regions only)
    B_hot = planck_lambda(T_hot, wavelength)

    # Wavelength-dependent transit depth
    delta_lambda = delta_gray * (B_hot / F_star)

    return delta_lambda
