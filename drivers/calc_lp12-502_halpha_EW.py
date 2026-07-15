"""
Analog of `calc_tic3006_halpha_EW.py`, but for LP 12-502 (TIC 402980664), using
the HIRES specriver cache from the night of 2024 Dec 03 (run label j546), in
`results/movie_sixpanel_specriver/TIC_402980664_j546_science_wob/`.

Only the first 3.5 hours of observations are used: a flare contaminates the
Halpha profile toward the end of the night (see the specriver river / mosaic in
`drivers/plot_hires_lp12-502_mosaic.py`).  The mean Halpha line profile is
recomputed from those first-3.5h spectra (rather than using the cached
`smooth_mean_flux`, which averages in the flare-contaminated tail).

As in the TIC 300651846 analysis, the mean profile is fit with gaussian,
lorentzian, pseudo-voigt, double gaussian, and double lorentzian models.  We
report reduced chi2, delta BIC, and best-fit parameters, then compute an
equivalent width by integrating the best (lowest-BIC) model.

Notable differences from `calc_tic3006_halpha_EW.py`:
  - Single night, so there is no multi-night loop / comparison plot.
  - The cache stores the velocity axis in km/s; here it is normalized by
    ADOPTED_VEQ = 24.2 km/s to reproduce the reference's "v_eq units" axis and
    its EW conversion factor.
  - The LP 12-502 line is narrow (confined to |dv/v_eq| < ~0.15) and the mean
    profile is unsmoothed, so the flux uncertainty is estimated from the
    continuum-wing scatter (std for |dv/v_eq| > 0.3) rather than from the
    smoothed-curve `p2p_rms x 50` heuristic used for TIC 300651846.
"""
from glob import glob
from os.path import join
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

from aesthetic.plot import set_style

from complexrotators.paths import RESULTSDIR

DATADIR = join(RESULTSDIR, 'movie_sixpanel_specriver')
OUTDIR = join(RESULTSDIR, 'LP12502_Halpha_static_EW')

# LP 12-502 / TIC 402980664 (per plot_hires_lp12-502_mosaic.py)
TICID = '402980664'
DATESTR = 'j546'
DATELABEL = '2024-12-03'
ADOPTED_VEQ = 24.2  # km/s, = 2*pi*Rstar/Prot
HALPHA_AIR_WAV = 6562.8  # Angstrom
C_KMS = 299792.458

# only use spectra taken within this many hours of the first observation; a
# flare contaminates the Halpha profile afterward.
MAX_ELAPSED_HR = 3.5

# velocity-domain (in v_eq units) selections
WING_ABS_DV = 0.3    # |dv| > this defines the continuum used for the noise est.
NARROW_ABS_DV = 0.15  # |dv| < this defines the "narrow" EW integration domain


def find_pklpath(datestr: str) -> str:
    """Locate the specriver cache pickle for the LP 12-502 night.

    Args:
        datestr: run label, e.g. 'j546'.

    Returns:
        str: path to the matching *specriver*cache*pkl.
    """
    pattern = join(DATADIR, f'TIC_{TICID}_{datestr}_science_wob',
                   '*specriver*cache*pkl')
    matches = glob(pattern)
    if len(matches) != 1:
        raise ValueError(f'Expected 1 match for {pattern}, got {matches}')
    return matches[0]


def load(pklpath: str) -> tuple:
    """Load the mean Halpha profile over the first `MAX_ELAPSED_HR` hours.

    The cache stores per-spectrum continuum-normalized profiles in
    `orig_yvals` (a list of 1D arrays), the shared velocity axis in `xvals`
    (km/s), and observation times in `spectimes` (BTJD).  We select the
    spectra within `MAX_ELAPSED_HR` of the first observation, average them,
    normalize the velocity axis by `ADOPTED_VEQ`, and estimate the flux
    uncertainty from the continuum-wing scatter.

    Args:
        pklpath: path to a specriver cache pickle.

    Returns:
        tuple: relative velocity in v_eq units (dv), mean normalized flux
            (flx), estimated per-pixel flux uncertainty (flx_err), and the
            number of spectra averaged (n_spec).
    """
    with open(pklpath, 'rb') as f:
        d = pickle.load(f)

    spectimes = np.asarray(d['spectimes'])
    orig_yvals = np.asarray(d['orig_yvals'])  # (n_spec, n_wav)
    xvals = np.asarray(d['xvals'])[0]         # km/s, shared across spectra

    elapsed_hr = (spectimes - spectimes.min()) * 24.0
    sel_time = elapsed_hr <= MAX_ELAPSED_HR
    n_spec = int(sel_time.sum())

    flx = np.nanmean(orig_yvals[sel_time], axis=0)
    dv = xvals / ADOPTED_VEQ

    # noise estimate: scatter of the mean profile in the line-free continuum.
    wing = np.abs(dv) > WING_ABS_DV
    flx_err = float(np.nanstd(flx[wing]))

    return dv, flx, flx_err, n_spec


def gauss_model(x: np.ndarray, amp: float, center: float, sigma: float,
                offset: float) -> np.ndarray:
    """Gaussian model for an emission line.

    Args:
        x (np.ndarray): Array of velocities.
        amp (float): Amplitude of the Gaussian.
        center (float): Center position of the Gaussian.
        sigma (float): Gaussian width (standard deviation).
        offset (float): Baseline offset.

    Returns:
        np.ndarray: Gaussian function evaluated at x.
    """
    return offset + amp * np.exp(-((x - center) ** 2) / (2 * sigma ** 2))


def lorentz_model(x: np.ndarray, amp: float, center: float, gamma: float,
                  offset: float) -> np.ndarray:
    """Lorentzian model for an emission line.

    Args:
        x (np.ndarray): Array of velocities.
        amp (float): Amplitude of the Lorentzian.
        center (float): Center position of the Lorentzian.
        gamma (float): Half-width at half-maximum.
        offset (float): Baseline offset.

    Returns:
        np.ndarray: Lorentzian function evaluated at x.
    """
    return offset + amp * (gamma ** 2 / ((x - center) ** 2 + gamma ** 2))


def pseudo_voigt(x: np.ndarray, amp: float, center: float, sigma: float,
                 eta: float, offset: float) -> np.ndarray:
    """Pseudo-Voigt profile: a linear combination of Gaussian and Lorentzian.

    Args:
        x (np.ndarray): Array of velocities.
        amp (float): Amplitude of the profile.
        center (float): Center position of the profile.
        sigma (float): Width parameter for both Gaussian and Lorentzian.
        eta (float): Mixing parameter (0: pure Gaussian, 1: pure Lorentzian).
        offset (float): Baseline offset.

    Returns:
        np.ndarray: Pseudo-Voigt function evaluated at x.
    """
    gauss = np.exp(-((x - center) ** 2) / (2 * sigma ** 2))
    lorentz = sigma ** 2 / ((x - center) ** 2 + sigma ** 2)
    return offset + amp * (eta * lorentz + (1 - eta) * gauss)


def double_gauss_model(x: np.ndarray, amp1: float, center1: float,
                       sigma1: float, amp2: float, center2: float,
                       sigma2: float, offset: float) -> np.ndarray:
    """Double Gaussian model for an emission line.

    Args:
        x (np.ndarray): Array of velocities.
        amp1 (float): Amplitude of the first Gaussian.
        center1 (float): Center of the first Gaussian.
        sigma1 (float): Width of the first Gaussian.
        amp2 (float): Amplitude of the second Gaussian.
        center2 (float): Center of the second Gaussian.
        sigma2 (float): Width of the second Gaussian.
        offset (float): Baseline offset.

    Returns:
        np.ndarray: Sum of two Gaussian functions with an offset.
    """
    gauss1 = amp1 * np.exp(-((x - center1) ** 2) / (2 * sigma1 ** 2))
    gauss2 = amp2 * np.exp(-((x - center2) ** 2) / (2 * sigma2 ** 2))
    return offset + gauss1 + gauss2


def double_lorentz_model(x: np.ndarray, amp1: float, center1: float,
                         gamma1: float, amp2: float, center2: float,
                         gamma2: float, offset: float) -> np.ndarray:
    """Double Lorentzian model for an emission line.

    Args:
        x (np.ndarray): Array of velocities.
        amp1 (float): Amplitude of the first Lorentzian.
        center1 (float): Center of the first Lorentzian.
        gamma1 (float): Half-width at half-maximum of the first Lorentzian.
        amp2 (float): Amplitude of the second Lorentzian.
        center2 (float): Center of the second Lorentzian.
        gamma2 (float): Half-width at half-maximum of the second Lorentzian.
        offset (float): Baseline offset.

    Returns:
        np.ndarray: Sum of two Lorentzian functions with an offset.
    """
    lorentz1 = amp1 * (gamma1 ** 2 / ((x - center1) ** 2 + gamma1 ** 2))
    lorentz2 = amp2 * (gamma2 ** 2 / ((x - center2) ** 2 + gamma2 ** 2))
    return offset + lorentz1 + lorentz2


def fit_model(model_func, dv: np.ndarray, flx: np.ndarray,
              flx_err: float, p0: list, bounds: tuple) -> tuple:
    """Fit a given model to the spectrum data.

    Args:
        model_func (callable): The model function to fit.
        dv (np.ndarray): Array of relative velocities.
        flx (np.ndarray): Array of normalized flux values.
        flx_err (float): Estimated uncertainty in flux.
        p0 (list): Initial guess for the model parameters.
        bounds (tuple): Lower and upper bounds for the parameters.

    Returns:
        tuple: Optimal parameters, covariance matrix, fitted flux, number
            of free parameters, reduced chi-squared, and BIC. If the fit
            fails, returns None for parameters and covariance, a NaN array
            for fitted flux, and np.nan for statistics.
    """
    sigma_arr = np.full_like(flx, flx_err)
    try:
        popt, pcov = curve_fit(
            model_func, dv, flx, p0=p0, bounds=bounds,
            sigma=sigma_arr, absolute_sigma=True, maxfev=10000)
        fit_flux = model_func(dv, *popt)
        resid = (flx - fit_flux) / sigma_arr
        chi2 = np.sum(resid ** 2)
        dof = len(dv) - len(popt)
        red_chi2 = chi2 / dof
        bic = chi2 + len(popt) * np.log(len(dv))
    except RuntimeError:
        popt, pcov = None, None
        fit_flux = np.full_like(flx, np.nan)
        red_chi2 = np.nan
        bic = np.nan
    n_params = len(popt) if popt is not None else np.nan
    return popt, pcov, fit_flux, n_params, red_chi2, bic


def analyze_night(pklpath: str, datestr: str, outdir: str) -> dict:
    """Fit the LP 12-502 mean Halpha emission line with various functional
    forms and compute the equivalent width.

    Loads the mean spectrum over the first `MAX_ELAPSED_HR` hours, fits it
    with Gaussian, Lorentzian, Pseudo-Voigt, Double Gaussian, and Double
    Lorentzian profiles, makes a two-panel fits/residuals plot, prints model
    comparison and best-fit parameter tables, then identifies the best model
    by BIC and integrates it to get the equivalent width.

    Args:
        pklpath: path to the night's specriver cache pickle.
        datestr: run label, e.g. 'j546', used for labeling outputs.
        outdir: directory to write the plot and result CSVs to.

    Returns:
        dict: summary of the best-fit model and EW measurements.
    """
    dv, flx, flx_err, n_spec = load(pklpath)
    print(f'Averaged {n_spec} spectra taken within {MAX_ELAPSED_HR} h of the '
          f'first observation.')

    amp_guess = np.max(flx) - np.min(flx)
    center_guess = 0
    width_guess = 0.1  # v_eq units; the line is narrow (|dv| < ~0.15)
    offset_guess = np.median(flx[np.abs(dv) > WING_ABS_DV])

    bounds_gauss = (
        [0, np.min(dv), 0, -np.inf],
        [np.inf, np.max(dv), np.ptp(dv), np.inf]
    )
    bounds_lorentz = bounds_gauss
    bounds_voigt = (
        [0, np.min(dv), 0, 0, -np.inf],
        [np.inf, np.max(dv), np.ptp(dv), 1, np.inf]
    )
    bounds_dgauss = (
        [0, np.min(dv), 0, 0, np.min(dv), 0, -np.inf],
        [np.inf, np.max(dv), np.ptp(dv), np.inf, np.max(dv),
         np.ptp(dv), np.inf]
    )
    bounds_dlorentz = bounds_dgauss

    models = {
        'Gaussian': {
            'func': gauss_model,
            'p0': [amp_guess, center_guess, width_guess, offset_guess],
            'bounds': bounds_gauss
        },
        'Lorentzian': {
            'func': lorentz_model,
            'p0': [amp_guess, center_guess, width_guess, offset_guess],
            'bounds': bounds_lorentz
        },
        'Pseudo-Voigt': {
            'func': pseudo_voigt,
            'p0': [amp_guess, center_guess, width_guess, 0.5, offset_guess],
            'bounds': bounds_voigt
        },
        'Double Gaussian': {
            'func': double_gauss_model,
            'p0': [amp_guess * 0.6, center_guess, width_guess,
                   amp_guess * 0.4, center_guess, width_guess * 3,
                   offset_guess],
            'bounds': bounds_dgauss
        },
        'Double Lorentzian': {
            'func': double_lorentz_model,
            'p0': [amp_guess * 0.6, center_guess, width_guess,
                   amp_guess * 0.4, center_guess, width_guess * 3,
                   offset_guess],
            'bounds': bounds_dlorentz
        }
    }

    fit_results = {}
    summary_list = []

    for name, info in models.items():
        popt, pcov, fit_flux, n_params, red_chi2, bic = fit_model(
            info['func'], dv, flx, flx_err, info['p0'], info['bounds'])
        fit_results[name] = {
            'popt': popt,
            'pcov': pcov,
            'fit_flux': fit_flux,
            'n_params': n_params,
            'red_chi2': red_chi2,
            'BIC': bic
        }
        summary_list.append({
            'Model': name,
            'n_params': n_params,
            'Reduced Chi2': red_chi2,
            'BIC': bic
        })

    summary_df = pd.DataFrame(summary_list)
    # delta BIC relative to the best (lowest-BIC) model.
    if summary_df['BIC'].notna().any():
        summary_df['delta_BIC'] = summary_df['BIC'] - summary_df['BIC'].min()

    # x-axis in wavelength: dv (v_eq units) -> Angstroms about H-alpha.
    conversion_factor = (ADOPTED_VEQ * HALPHA_AIR_WAV) / C_KMS
    wav = HALPHA_AIR_WAV + dv * conversion_factor

    set_style('science')
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)
    ax1.errorbar(wav, flx, yerr=np.full_like(flx, flx_err),
                 fmt='k.', label='Data')
    for name, result in fit_results.items():
        ax1.plot(wav, result['fit_flux'], label=name)
    ax1.set_ylabel('Normalized Flux')
    ax1.legend()
    ax1.set_title(
        f'LP 12-502 ({datestr}, {DATELABEL}): Fits to H-alpha Emission Line')

    for name, result in fit_results.items():
        ax2.plot(wav, flx - result['fit_flux'], label=name)
    ax2.set_xlabel(r'Wavelength [$\mathrm{\AA}$]')
    ax2.set_ylabel('Residuals')
    ax2.legend()
    ax2.set_title('Residuals of the Fits')

    for ax in (ax1, ax2):
        ax.minorticks_on()
        ax.tick_params(which='both', top=True, right=True, direction='in')
        ax.tick_params(which='major', length=5)
        ax.tick_params(which='minor', length=3)

    plt.tight_layout()
    os.makedirs(outdir, exist_ok=True)
    figpath = join(outdir, f'LP12502_{datestr}_emission_line_fits.png')
    fig.savefig(figpath, dpi=300)
    plt.close(fig)
    print(f'Made {figpath}')

    print(f'\n=== {datestr} ({DATELABEL}) ===')
    print("Model Fit Statistics:")
    print(summary_df)

    # Build summary table for best-fit parameters and their uncertainties.
    param_names = {
        'Gaussian': ['amp', 'center', 'sigma', 'offset'],
        'Lorentzian': ['amp', 'center', 'gamma', 'offset'],
        'Pseudo-Voigt': ['amp', 'center', 'sigma', 'eta', 'offset'],
        'Double Gaussian': ['amp1', 'center1', 'sigma1',
                            'amp2', 'center2', 'sigma2', 'offset'],
        'Double Lorentzian': ['amp1', 'center1', 'gamma1',
                              'amp2', 'center2', 'gamma2', 'offset']
    }
    params_summary_list = []
    for model_name, result in fit_results.items():
        popt = result['popt']
        pcov = result['pcov']
        if popt is not None and pcov is not None:
            perr = np.sqrt(np.diag(pcov))
            for pname, val, err in zip(param_names[model_name],
                                       popt, perr):
                params_summary_list.append({
                    'Model': model_name,
                    'Parameter': pname,
                    'Value': val,
                    'Uncertainty': err
                })
    params_summary_df = pd.DataFrame(params_summary_list)
    print("\nBest-fit Parameters Summary:")
    print(params_summary_df)

    statspath = join(outdir, f'LP12502_{datestr}_model_fit_statistics.csv')
    summary_df.to_csv(statspath, index=False)
    parampath = join(outdir, f'LP12502_{datestr}_bestfit_parameters.csv')
    params_summary_df.to_csv(parampath, index=False)
    print(f'Made {statspath}\nMade {parampath}')

    # Identify the best model by lowest BIC.
    valid_models = {m: info['BIC'] for m, info in fit_results.items()
                    if not np.isnan(info['BIC'])}

    result_row = {
        'datestr': datestr,
        'datelabel': DATELABEL,
        'n_spec': n_spec,
        'flx_err': flx_err,
        'best_model': np.nan,
        'red_chi2': np.nan,
        'BIC': np.nan,
        'EW_narrow_A': np.nan,
        'EW_total_A': np.nan,
    }

    if valid_models:
        best_model_name = min(valid_models, key=valid_models.get)
        best_model_func = models[best_model_name]['func']
        best_params = fit_results[best_model_name]['popt']
        best_fit_flux = best_model_func(dv, *best_params)

        # conversion_factor (dv in v_eq units -> Angstroms) defined above.
        # Equivalent widths. Negative EW => net emission (flux > continuum).
        sel_domain = np.abs(dv) < NARROW_ABS_DV
        ew_narrow = np.trapezoid((best_fit_flux[sel_domain] - 1) *
                                 conversion_factor, dv[sel_domain])
        ew_total = np.trapezoid((best_fit_flux - 1) * conversion_factor, dv)

        ew_table = pd.DataFrame({
            'Domain': [f'|dv|<{NARROW_ABS_DV}', 'All dv'],
            'Equivalent Width (Angstroms)': [ew_narrow, ew_total]
        })

        print("\nEquivalent Width Calculations:")
        print("Best model (lowest BIC):", best_model_name)
        print(ew_table)

        result_row.update({
            'best_model': best_model_name,
            'red_chi2': fit_results[best_model_name]['red_chi2'],
            'BIC': fit_results[best_model_name]['BIC'],
            'EW_narrow_A': ew_narrow,
            'EW_total_A': ew_total,
        })
    else:
        print("No valid model fits to compute equivalent widths.")

    return result_row


def calc_lp12502_halpha_EW():
    """Run the H-alpha EW analysis for the LP 12-502 j546 night (2024 Dec 03).

    Uses only the first `MAX_ELAPSED_HR` hours of spectra to avoid a flare that
    contaminates the profile at the end of the night, fits the mean Halpha
    line, and writes the best-fit model and equivalent width to a summary CSV.
    """
    os.makedirs(OUTDIR, exist_ok=True)

    pklpath = find_pklpath(DATESTR)
    row = analyze_night(pklpath, DATESTR, OUTDIR)

    result_df = pd.DataFrame([row])
    summarypath = join(OUTDIR, 'LP12502_Halpha_EW_summary.csv')
    result_df.to_csv(summarypath, index=False)
    print(f'\nMade {summarypath}')
    print(result_df)


if __name__ == "__main__":
    calc_lp12502_halpha_EW()
