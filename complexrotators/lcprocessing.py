"""
Contents:
    | cpv_periodsearch
    | count_phased_local_minima
    | prepare_cpv_light_curve: retrieve all relevant data from SPOC/FITS LC
    | p2p_rms
"""

#######################################
# ASTROBASE IMPORTS CAN BREAK LOGGING #
#######################################
from astrobase.periodbase import pgen_lsp, stellingwerf_pdm
from astrobase.checkplot import checkplot_png
from astrobase.lcmath import phase_magseries, phase_bin_magseries

#############
## LOGGING ##
#############
import logging
from complexrotators import log_sub, log_fmt, log_date_fmt

DEBUG = False
if DEBUG:
    level = logging.DEBUG
else:
    level = logging.INFO
LOGGER = logging.getLogger(__name__)
logging.basicConfig(
    level=level,
    style=log_sub,
    format=log_fmt,
    datefmt=log_date_fmt,
    force=True
)

LOGDEBUG = LOGGER.debug
LOGINFO = LOGGER.info
LOGWARNING = LOGGER.warning
LOGERROR = LOGGER.error
LOGEXCEPTION = LOGGER.exception

#############
## IMPORTS ##
#############
import numpy as np, pandas as pd
from numpy import array as nparr

import os, multiprocessing, pickle
from complexrotators.paths import RESULTSDIR, DATADIR

from astropy.io import fits

from scipy.signal import find_peaks

from wotan import flatten
from copy import deepcopy

nworkers = multiprocessing.cpu_count()

def cpv_periodsearch(times, fluxs, starid, outdir, t0=None,
                     periodogram_method="pdm", runperiodsearch=1):
    """
    Given time and flux, run a period-search for objects expected to be complex
    rotators.

    A few plots and pickle files will be written to `outdir` using the `starid`
    string.

    Args:

        times (np.ndarray):
            Array of times.

        fluxs (np.ndarray):
            Array of fluxes.

        starid (str):
            Identifier used for cacheing.

        outdir (str):
            Path used for cacheing.

        t0 (None, str, int, or float):
            Epoch at which to phase.  None defaults to t0=1618.  Giving the
            string "binmin" defaults to phase-folding, and taking the
            arg-minimum.  Any int or float will be passed as the manual phase.

        periodogram_method (str):
            "pdm" (phase dispersion minimization) or "ls" (lomb-scargle).

    Returns:

        dict : results

            A dictionary of the results, containing:
                'lsp':periodogram results, 'fine_lsp':fine periodogram results,
                'times':times, 'fluxs':fluxs, 'period':fine_lsp['bestperiod'],
                't0': t0, 'outdir':outdir, 'periodogram_method': ...
            Note that just because the keys are "lsp", the actual method being
            used depends on periodogram_method
    """

    assert isinstance(starid, str)

    pklpath = os.path.join(outdir, f"{starid}_cpv_periodsearch.pkl")
    if os.path.exists(pklpath):
        LOGINFO(f"Found {pklpath}, loading and continuing.")
        with open(pklpath, 'rb') as f:
            d = pickle.load(f)
        return d

    sep = 1
    # eg., a single TESS sector at 2-minute cadence has 2e4 points.  this cuts
    # it down for period-search purposes to 4e3, which helps the runtime!
    if len(times) > 1e4:
        sep = 5
    # eg., a single TESS sector at 2-minute cadence has 1.2e5 points.
    if len(times) > 1e5:
        sep = 50

    startp, endp = 0.05, 10

    # for the fine-tuning
    delta_P = 0.2
    stepsize = 1e-5

    LOGINFO(f'Beginning period search for {starid}')

    if periodogram_method == 'ls' and runperiodsearch:

        lsp = pgen_lsp(
            times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
            startp=startp, endp=endp, autofreq=True, sigclip=5.0, nbestpeaks=10
        )

        fine_lsp = pgen_lsp(
            times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
            startp=lsp['bestperiod']-delta_P*lsp['bestperiod'],
            endp=lsp['bestperiod']+delta_P*lsp['bestperiod'],
            autofreq=False, sigclip=5.0, stepsize=stepsize
        )

    elif periodogram_method == 'pdm' and runperiodsearch:

        lsp = stellingwerf_pdm(
            times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
            startp=startp, endp=endp, autofreq=True, sigclip=5.0, nbestpeaks=10
        )

        if lsp['bestperiod'] < 2:
            fine_lsp = stellingwerf_pdm(
                times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
                startp=lsp['bestperiod']-delta_P*lsp['bestperiod'],
                endp=lsp['bestperiod']+delta_P*lsp['bestperiod'],
                autofreq=False, sigclip=5.0, stepsize=stepsize
            )
        else:
            LOGINFO(
                f"Found P={lsp['bestperiod']:.3f} d; skipping fine period "
                f"search."
            )
            fine_lsp = deepcopy(lsp)

    LOGINFO(42*'.')
    LOGINFO(f"Standard autofreq period: {lsp['bestperiod']:.7f} d")
    LOGINFO(f"Fine period: {fine_lsp['bestperiod']:.7f} d")
    LOGINFO(f"Fine - standard: {fine_lsp['bestperiod']-lsp['bestperiod']:.7f} d")
    LOGINFO(42*'.')

    outfile = os.path.join(
        outdir, f'{starid}_{periodogram_method}_subset_checkplot.png'
    )
    if not os.path.exists(outfile):
        try:
            checkplot_png(lsp, times, fluxs, fluxs*1e-4, magsarefluxes=True,
                          phasewrap=True, phasesort=True, phasebin=0.002,
                          minbinelems=7, plotxlim=(-0.6,0.6), plotdpi=75,
                          outfile=outfile, verbose=True)
        except Exception as e:
            LOGEXCEPTION(e)
            LOGINFO("ontinuing...")
            pass

    if t0 is None:
        # default phase
        t0 = 1642.

    elif t0 == 'binmin':
        # bin the phase-fold to 50 points, take the minimum index.

        period = fine_lsp['bestperiod']
        x,y = times, fluxs-np.nanmean(fluxs)
        t0_ini = np.nanmin(x)
        _pd = phase_magseries(x, y, period, t0_ini, wrap=False,
                              sort=False)
        x_fold = _pd['phase']
        y = _pd['mags']
        bs_days = period/50
        try:
            orb_bd = phase_bin_magseries(x_fold, y, binsize=bs_days, minbinelems=3)
            min_phase = orb_bd['binnedphases'][np.argmin(orb_bd['binnedmags'])]
            t0 = t0_ini + min_phase*period
        except (ValueError, TypeError):
            # can be raised for very short periods...
            t0 = t0_ini

    elif isinstance(t0, (int, float)):
        pass

    else:
        raise NotImplementedError

    d = {
        'lsp':lsp, 'fine_lsp':fine_lsp, 'times':times, 'fluxs':fluxs,
        'period':fine_lsp['bestperiod'], 't0': t0, 'outdir':outdir,
        'periodogram_method': periodogram_method
        }

    with open(pklpath, 'wb') as f:
        pickle.dump(d, f)
        LOGINFO(f'Made {pklpath}')

    return d


def count_phased_local_minima(
    time, flux, t0, period,
    method="medianfilt_findpeaks",
    binsize_phase_units=0.02,
    height='2_P2P',
    width=2,
    window_length_phase_units=0.1,
    max_splines=None,
    height_limit=1e-3,
    pre_normalize=False
    ):
    """
    Given time, flux, epoch, and period, phase the light curve and count the
    number of "dips" (local minima) in the phase fold.

    Args:
        method (str):
            One of ["findpeaks", "medianfilt_findpeaks",
            "psplinefilt_findpeaks", "sinefilt_findpeaks"].  These are
            discussed further below.

        binsize_phase_units (float):
            Size of bins in units of phase.  e.g., 0.01 corresponds to 100
            points.

        height (float or str):
            Minimum height of peak required for it to be signficant. 1e-3 means
            0.1% in flux.  "1_MAD" means 1*median absolute deviation.  This
            needs to be a "_" separated string.  "2_P2P" means 2*the
            point-to-point rms.
            Actual adopted values will be max([height, height_limit]), so in
            practice dips smaller than height_limit will never be acquired.

        width (float):
            Minimum width of peak in units of samples.  E.g., 3 with
            binsize_phase_units of 0.01 means "at least 0.03 in phase".

        window_length_phase_units (float):
            Used only if method includes a windowed-slider, in which case this
            will be the window length in units of phase.

        max_splines (None or int):
            Needed if method == "psplinefilt_findpeaks".

        height_limit (float):
            See above.

    Returns:
        dictionary of results includes the number of peaks, their widths, and
        their heights.

    NOTES:
    Both the "findpeaks" and "medianfilt_findpeaks" methods rely on
    scipy.signal.find_peaks (link below).  The latter incorporates an initial
    step of median-smoothing, which is useful for removing smooth (presumably
    spot-induced) variability.

    LINKS:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html#scipy.signal.find_peaks
    """

    x,y = time, flux
    _pd = phase_magseries(x, y, period, t0, wrap=False, sort=True)

    DO_WRAP = 1
    if DO_WRAP:
        _pd['phase'] = np.concatenate(
            (_pd['phase']-1.0, _pd['phase'], _pd['phase']+1.0)
        )
        _pd['mags'] = np.concatenate(
            (_pd['mags'], _pd['mags'], _pd['mags'])
        )
        if isinstance(max_splines, int):
            max_splines *= 3

    x_fold = _pd['phase']
    y = _pd['mags']

    orb_bd = phase_bin_magseries(
        x_fold, y, binsize=binsize_phase_units, minbinelems=3
    )

    # x and flat_y are the points to search.  these are already phase-ordered.
    x, _y = orb_bd['binnedphases'], orb_bd['binnedmags']
    p2p_raw = p2p_rms(_y)
    a_95_5 = np.nanpercentile(_y, 95) - np.nanpercentile(_y, 5)
    a_max_min = np.nanmax(_y) - np.nanmin(_y)
    _y_mean = np.nanmean(_y)

    if pre_normalize:
        # normalize to a 1% semi-amplitude signal.
        a_desired = 0.02
        mult_fac = a_max_min / a_desired
        # rescale it
        _y_preflat = _y_mean + (_y - _y_mean)/mult_fac
    else:
        _y_preflat = _y
        mult_fac = 1


    nsplines = None

    if method == 'findpeaks':
        trend_y = None
        flat_y = _y_preflat * 1.

    elif method == 'medianfilt_findpeaks':
        flat_y, trend_y = flatten(
            x, _y_preflat, method='median', return_trend=True,
            break_tolerance=1, window_length=window_length_phase_units,
            edge_cutoff=1e-3
        )

    elif method == 'sinefilt_findpeaks':
        flat_y, trend_y = flatten(
            x, _y_preflat, method='cosine', robust=True, return_trend=True,
            break_tolerance=1, window_length=window_length_phase_units
        )

    elif method == 'psplinefilt_findpeaks':
        flat_y, trend_y, nsplines = flatten(
            x, _y_preflat, method='pspline', max_splines=max_splines,
            edge_cutoff=1e-3, stdev_cut=2, return_trend=True,
            return_nsplines=True, verbose=True
        )

        # estimate p2p_rms on a heavily whitened light curve:
        # 33*3 splines corresponds to a window per phase of φ=0.033 -> heavily
        # whitened!
        veryflat_y, _ = flatten(
            x, _y_preflat, method='pspline', max_splines=33*3,
            edge_cutoff=1e-3, stdev_cut=5, return_trend=True,
            return_nsplines=False, verbose=True
        )

    if method == 'psplinefilt_findpeaks':
        p2p_est = p2p_rms(veryflat_y)
    else:
        p2p_est = p2p_raw * 1.


    offset = np.nanmean(flat_y)
    flat_y -= offset

    # number of points in one full cycle
    N = int(1/binsize_phase_units)

    mad = np.nanmedian( np.abs( flat_y - np.nanmedian(flat_y)  ) )

    if isinstance(height, float):
        height = height
    if isinstance(height, str):
        if 'MAD' in height:
            height = mad * float(height.split("_")[0])
        if 'P2P' in height:
            height = p2p_est * float(height.split("_")[0])

    # say the p2p estimate is very small.  in such cases, we do not really care
    # about the dips... 
    height = max([ height, height_limit ])

    LOGINFO(f"p2p_raw {p2p_raw:.1e}, p2p_est {p2p_est:.1e}, mad {mad:.1e}, height {height:.1e}")
    peaks, properties = find_peaks(
        -flat_y, height=height, width=width, rel_height=0.5
    )

    #import matplotlib.pyplot as plt
    #plt.scatter(x, flat_y)
    #plt.savefig('temp.png')
    #import IPython; IPython.embed()
    #assert 0

    # construct the list of unique peaks, since there are duplicates because of
    # the need to phase wrap.
    #
    # For example,
    #
    # unique_peak_indices = {2, 67, 76, 84}
    # peaks = array([ 67,  84, 102, 167, 176, 184, 202, 267, 284])
    # 
    # will yield
    #
    # peak_list = [67, 84, 2, 76]
    # sel = [T, T, T, F, T, F, F, F, F]

    unique_peak_indices = np.array(sorted(list(set(peaks % N))))

    # first, remove any peaks within 1 or 2 phase units of each other, and keep the
    # second of each such pair
    # e.g., peaks = [50, 100, 150, 201]
    # yields unique_peak_indices = [0, 1, 50].
    # here, we would drop the "0" case.
    if len(unique_peak_indices) > 1:
        _sel = np.ones(len(unique_peak_indices)).astype(bool)
        if not np.any(np.diff(unique_peak_indices) == 1):
            pass
        else:
            _sel = np.roll(unique_peak_indices, -1) - unique_peak_indices != 1
            _sel &= np.roll(unique_peak_indices, -1) - unique_peak_indices != 2
        unique_peak_indices = unique_peak_indices[_sel]

    sel = np.zeros(len(peaks)).astype(bool)
    peak_list = []
    for ix, peak in enumerate(peaks):
        if peak % N in unique_peak_indices and peak % N not in peak_list:
            sel[ix] = True
            peak_list.append(peak % N)
    assert np.sum(sel) == len(unique_peak_indices)

    peaks = peaks[sel]
    for k, v in properties.items():
        properties[k] = v[sel]
    params = "left_bases,right_bases,widths,left_ips,right_ips".split(",")
    for param in params:
        properties[param+"_phaseunits"] = properties[param]*binsize_phase_units

    if not pre_normalize:
        binned_trend_flux = trend_y
    else:
        binned_trend_flux = (
            mult_fac * (trend_y - _y_mean) + _y_mean
        )

    r = {
        'N_peaks': len(peaks),
        'peaks_phaseunits': np.mod(peaks*binsize_phase_units, 1),
        'peaks': peaks,
        'properties': properties,
        'time': time,
        'flux': flux,
        't0': t0,
        'period': period,
        'binsize_phase_units': binsize_phase_units,
        'height': height,
        'width': width,
        'phase': _pd['phase'],
        'phase_flux': _pd['mags'],
        'binned_phase': x,
        'binned_search_flux': flat_y,
        'p2p_raw': p2p_raw,
        'p2p_est': p2p_est,
        'a_95_5': a_95_5,
        'mad': mad,
        'binned_orig_flux': _y,
        # if pre_normalize is True, then _y_preflat is what is used to actually
        # fit the penalized spline and detrend.
        # _y_preflat = _y_mean + (_y - _y_mean)/mult_fac
        'binned_preflat_flux': _y_preflat,
        'mult_fac': mult_fac,
        '_y_mean': _y_mean,
        'binned_trend_flux': binned_trend_flux,
        'nsplines_total': nsplines,
        'nsplines_singlephase': nsplines/3
    }

    return r


def prepare_cpv_light_curve(lcpath, cachedir, returncadenceno=0,
                            lcpipeline='spoc2min'):
    """
    Given a light curve (SPOC 2-minute or QLP), remove non-zero quality flags,
    median-normalize, and run a 5-day median filter over the light curve.
    Cache the output.
    """

    hl = fits.open(lcpath)
    hdr = hl[0].header
    d = hl[1].data

    # metadata
    sector = hdr["SECTOR"]
    ticid = hdr["TICID"]

    # light curve data
    TIMEKEYDICT = {
        'qlp': 'TIME',
        'spoc2min': 'TIME',
        'cdips': 'TMID_BJD'
    }
    time = d[TIMEKEYDICT[lcpipeline]]
    FLUXKEYDICT = {
        'spoc2min': 'PDCSAP_FLUX',
        # As of 11/21/2023, the QLP has switched from "KSPSAP_FLUX" to
        # "DET_FLUX".  This is because they "changed their detrending
        # algorithm".  Qualitatively similar flattening in the latter.
        'qlp': ['KSPSAP_FLUX', 'DET_FLUX'],
        'cdips': 'PCA3'
    }
    if lcpipeline in ['spoc2min', 'cdips']:
        flux = d[FLUXKEYDICT[lcpipeline]]
    elif lcpipeline == 'qlp':
        if 'KSPSAP_FLUX' in d.names:
            flux = d['KSPSAP_FLUX']
        elif 'DET_FLUX' in d.names:
            flux = d['DET_FLUX']
        else:
            raise NotImplementedError

    if lcpipeline == 'cdips':
        from cdips.utils.lcutils import _given_mag_get_flux
        flux = _given_mag_get_flux(flux)

    QUALITYKEYDICT = {
        'spoc2min': 'QUALITY',
        'qlp': 'QUALITY',
        'cdips': 'IRQ3'
    }

    qual = d[QUALITYKEYDICT[lcpipeline]]

    if lcpipeline in ['spoc2min', 'qlp']:
        cadenceno = d['CADENCENO']
        sel = (qual == 0)

    elif lcpipeline == 'cdips':
        sel = (qual == 'G')

    # remove non-zero quality flags
    x_obs = time[sel]
    y_obs = flux[sel]
    if lcpipeline in ['spoc2min', 'qlp']:
        cadenceno_obs = cadenceno[sel]

    if np.isfinite(y_obs).sum() < 20:
        return (None, None, None, None, None, None, None, None,
                None, None, None)

    # normalize around 1
    y_obs /= np.nanmedian(y_obs)

    # NOTE: you could consider removing flares using a time-windowed slider
    # here.  however, for purposes of finding the periods, they are a small
    # enough fraction of the duty cycle that they can probably be ignored.

    # what is the cadence?
    cadence_sec = int(np.round(np.nanmedian(np.diff(x_obs))*24*60*60))

    starid = f'{ticid}_S{str(sector).zfill(4)}_{cadence_sec}sec_{lcpipeline}'

    #
    # "light" detrending by default. (& cache it)
    #
    pklpath = os.path.join(cachedir, f"{starid}_dtr_lightcurve.pkl")
    if os.path.exists(pklpath):
        LOGINFO(f"Found {pklpath}, loading and continuing.")
        with open(pklpath, 'rb') as f:
            lcd = pickle.load(f)
        y_flat = lcd['y_flat']
        y_trend = lcd['y_trend']
        x_trend = lcd['x_trend']
    else:
        y_flat, y_trend = flatten(x_obs, y_obs, window_length=5.0,
                                  return_trend=True, method='median')
        x_trend = deepcopy(x_obs)
        lcd = {'y_flat':y_flat, 'y_trend':y_trend, 'x_trend':x_trend }
        with open(pklpath, 'wb') as f:
            pickle.dump(lcd, f)
            LOGINFO(f'Made {pklpath}')

    assert len(y_obs) == len(y_flat) == len(x_obs)

    if returncadenceno:
        assert lcpipeline in ['spoc2min', 'qlp']
        return (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend,
                cadence_sec, sector, starid, cadenceno_obs)
    else:
        return (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend,
                cadence_sec, sector, starid)


def p2p_rms(flux):
    """
    Calculate the 68th percentile of the distribution of the residuals from the
    median value of δF_i = F_{i} - F_{i+1}, where i is an index over time.
    """
    dflux = np.diff(flux)
    med_dflux = np.nanmedian(dflux)

    up_p2p = (
        np.nanpercentile( dflux-med_dflux, 84 )
        -
        np.nanpercentile( dflux-med_dflux, 50 )
    )
    lo_p2p = (
        np.nanpercentile( dflux-med_dflux, 50 )
        -
        np.nanpercentile( dflux-med_dflux, 16 )
    )

    p2p = np.nanmean([up_p2p, lo_p2p])

    return p2p

