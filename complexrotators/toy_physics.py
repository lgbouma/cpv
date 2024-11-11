"""
Contents:

Assuming dust cloud:
    dust_transit_depth
    dust_transit_depth_detailedopacity

Assuming plasma cloud + transit light source effect:
    planck_lambda
    stellar_flux
    transit_depth_tlse
"""
import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy.constants import h, c, k

from gasvsdust.synthetic_spectra import get_btsettl_spectrum

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


def dust_transit_depth_detailedopacity(
    wavelength: np.ndarray,
    r_star: float,
    Sigma: float,
    A_clump_ratio: float,
    opacitysource = 'Draine2003',
    opacitysubtype = None,
    R_V = 3.1
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
        opacitysource (str):
            "Draine2003",  ISM from Draine's website
            "G23",  Gordon2023 ISM
            "G29_pl",  Powerlaw+ G29-38 Alycia Weinberger fit
            "G29_cr",  Powerlaw+ G29-38 Alycia Weinberger with cranked up crystalline olivine
            "DL03", "C11" or "J13".  (Per gasvsdust.getters)
        opacitysubtype (str): Pick a specific dust composition from the DustEM suite:
            C11: PAH0_MC10, PAH1_MC10, amCBEx, amCBEx_2, aSilx
            DL01: PAH0_DL01, PAH1_DL01, Gra, Gra_2, aSil
            J13: CM20, CM20_2, aPyM5, aOlM5
            e.g. if you wanted the C11 silicate, set as "aSilx".

    Returns:
        np.ndarray: Dip depth at each wavelength.
    """
    # Star area
    A_star = np.pi * (r_star * 6.96e8)**2  # Convert solar radii to meters

    assert opacitysource in [
        "Draine2003", "DL03", "C11", "J13", "G23", "G29_pl", "G29_cr"
    ]

    from scipy.interpolate import interp1d
    if opacitysource == 'Draine2003':
        from gasvsdust.getters import get_draine_2003_astronsilicate_data
        df = get_draine_2003_astronsilicate_data(R_V=R_V)
        df = df.sort_values(by='lambda')
        xkey, ykey = 'lambda', 'K_abs'
        # Absorption + scattering  cross section per mass of dust [opacity] in cm^2/g
        fn = interp1d(df[xkey], df[ykey])
        kappa = fn(wavelength)

    elif opacitysource == 'G23':
        from dust_extinction.parameter_averages import G23
        from astropy import units as u
        extinction_model = G23(Rv=R_V)
        # Already interpolates
        A_lambda_over_Av = extinction_model(wavelength*u.micron)
        # Normalize to silicate 10um feature at R_V=3.1
        fudge = 3.65e4
        # Absorption + scattering cross section per mass of dust [opacity] in cm^2/g
        kappa = fudge*A_lambda_over_Av

    elif opacitysource in ['G29_pl', 'G29_cr']:
        from gasvsdust.getters import get_alycia_model_plus_powerlaw
        if opacitysource == 'G29_pl':
            model_id = 'g29_38_like_1200K_model'
        elif opacitysource == 'G29_cr':
            model_id = 'g29_38_like_more_crystalline'
        kappa = get_alycia_model_plus_powerlaw(wavelength, model_id)

    else:
        from gasvsdust.getters import get_DustEM_dat
        df = get_DustEM_dat(ref=opacitysource)
        df = df.sort_values(by='lambda')
        xkey, ykey = 'lambda', 'ext_tot'
        if opacitysubtype is not None:
            ykey = f'abs_{opacitysubtype}'
        # Absorption cross section per mass of dust [opacity] in cm^2/g
        # NOTE: should add scattering...
        fn = interp1d(df[xkey], df[ykey])
        kappa = fn(wavelength)

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
    B_lambda = (
        (2.0 * h * c**2) / (wavelength_m**5)
        /
        (np.exp((h * c) / (wavelength_m * k * T)) - 1)
    )

    return B_lambda


def stellar_flux(
    f: float,
    T_hot: float,
    T_cool: float,
    wavelength: np.ndarray,
    starmodel: str = 'blackbody'
    ) -> np.ndarray:
    """
    Calculate the total stellar flux accounting for cool regions and hot regions.

    Args:
        f (float): Fraction of the stellar surface covered by cool regions.
        T_hot (float): Temperature of the hot regions in Kelvin.
        T_cool (float): Temperature of the cool regions in Kelvin.
        wavelength (np.ndarray): Wavelength in microns.
        starmodel (str): 'blackbody' or 'btsettl' (BT-Settl CIFIST)

    Returns:
        np.ndarray: Total stellar flux at each wavelength.
    """
    if starmodel == 'blackbody':
        get_fn = planck_lambda
        use_wavelength = 1.* wavelength
    elif starmodel == 'btsettl':
        get_fn = get_btsettl_spectrum
        use_wavelength = 1e4 * wavelength

    B_hot = get_fn(T_hot, use_wavelength)
    B_cool = get_fn(T_cool, use_wavelength)
    #if starmodel == 'btsettl':
    #    import IPython; IPython.embed()

    # Total flux with fractional contributions from cool regions and hot regions
    F_star = (1 - f) * B_hot + f * B_cool

    return F_star


def transit_depth_tlse(
    f: float,
    T_hot: float,
    T_cool: float,
    delta_gray: float,
    wavelength: np.ndarray,
    starmodel: str = 'blackbody'
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
        starmodel (str): 'blackbody' or 'btsettl' (BT-Settl CIFIST)

    Returns:
        np.ndarray: Dip depth at each wavelength.
    """
    assert starmodel in ['blackbody', 'btsettl']

    # Stellar flux as seen from the entire star (hot regions + cool regions)
    F_star = stellar_flux(f, T_hot, T_cool, wavelength, starmodel=starmodel)

    # Flux from the transit chord (hot regions only)
    if starmodel == 'blackbody':
        B_hot = planck_lambda(T_hot, wavelength)
    elif starmodel == 'btsettl':
        B_hot = get_btsettl_spectrum(T_hot, wavelength*1e4)

    # Wavelength-dependent transit depth
    delta_lambda = delta_gray * (B_hot / F_star)

    return delta_lambda
