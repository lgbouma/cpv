"""
Contents:
    | cpv_periodsearch
    | count_phased_local_minima
    | prepare_cpv_light_curve
"""
import numpy as np, pandas as pd
from numpy import array as nparr

import os, multiprocessing, pickle
from complexrotators.paths import RESULTSDIR, DATADIR

from astropy.io import fits
from astrobase import periodbase, checkplot

from astrobase.lcmath import phase_magseries, phase_bin_magseries

from scipy.signal import find_peaks

from wotan import flatten
from copy import deepcopy


nworkers = multiprocessing.cpu_count()

def cpv_periodsearch(times, fluxs, starid, outdir, t0=None,
                     periodogram_method="pdm"):
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
        print(f"Found {pklpath}, loading and continuing.")
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

    startp, endp = 0.05, 5

    # for the fine-tuning
    delta_P = 0.2
    stepsize = 1e-5

    print(f'Beginning period search for {starid}')

    if periodogram_method == 'ls':

        lsp = periodbase.pgen_lsp(
            times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
            startp=startp, endp=endp, autofreq=True, sigclip=5.0
        )

        fine_lsp = periodbase.pgen_lsp(
            times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
            startp=lsp['bestperiod']-delta_P*lsp['bestperiod'],
            endp=lsp['bestperiod']+delta_P*lsp['bestperiod'],
            autofreq=False, sigclip=5.0, stepsize=stepsize
        )

    elif periodogram_method == 'pdm':

        lsp = periodbase.stellingwerf_pdm(
            times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
            startp=startp, endp=endp, autofreq=True, sigclip=5.0
        )

        fine_lsp = periodbase.stellingwerf_pdm(
            times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
            startp=lsp['bestperiod']-delta_P*lsp['bestperiod'],
            endp=lsp['bestperiod']+delta_P*lsp['bestperiod'],
            autofreq=False, sigclip=5.0, stepsize=stepsize
        )

    print(42*'.')
    print(f"Standard autofreq period: {lsp['bestperiod']:.7f} d")
    print(f"Fine period: {fine_lsp['bestperiod']:.7f} d")
    print(f"Fine - standard: {fine_lsp['bestperiod']-lsp['bestperiod']:.7f} d")
    print(42*'.')

    outfile = os.path.join(
        outdir, f'{starid}_{periodogram_method}_subset_checkplot.png'
    )
    checkplot.checkplot_png(lsp, times, fluxs, fluxs*1e-4,
                            magsarefluxes=True, phasewrap=True,
                            phasesort=True, phasebin=0.002, minbinelems=7,
                            plotxlim=(-0.6,0.6), plotdpi=200,
                            outfile=outfile, verbose=True)

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
        orb_bd = phase_bin_magseries(x_fold, y, binsize=bs_days, minbinelems=3)
        min_phase = orb_bd['binnedphases'][np.argmin(orb_bd['binnedmags'])]
        t0 = t0_ini + min_phase*period

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
        print(f'Made {pklpath}')

    return d


def count_phased_local_minima(
    time, flux, t0, period,
    binsize_phase_units=0.005,
    prominence=1e-3,
    width=6
    ):
    """
    Given time, flux, epoch, and period, phase the light curve and count the
    number of "dips" (local minima) in the phase fold.  Accomplish this as
    follows.

    First, bin to say 100 points per cycle.  Then invert the signal (multiply
    by minus one), and use scipy.signal.find_peaks.  Require a particular width
    in addition to the prominence for a signal to be called a "peak".  eg.,
    "width must be at least 0.03 in phase, amplitude must be at least 0.1% in
    flux", or similar.

    See bottom example at
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html#scipy.signal.find_peaks

    Args:
        binsize_phase_units (float):
            Size of bins in units of phase.  e.g., 0.01 corresponds to 100
            points.

        prominence (float):
            Minimum height of peak required for it to be signficant. 1e-3 means
            0.1% in flux.

        width (float):
            Minimum width of peak in units of samples.  E.g., 3 with
            binsize_phase_units of 0.01 means "at least 0.03 in phase".

    Returns:
        dictionary of results includes the number of peaks, their widths, and
        their prominences.
    """

    x,y = time, flux-np.nanmean(flux)
    _pd = phase_magseries(x, y, period, t0, wrap=True, sort=False)
    x_fold = _pd['phase']
    y = _pd['mags']

    orb_bd = phase_bin_magseries(
        x_fold, y, binsize=binsize_phase_units, minbinelems=5
    )
    min_phase = orb_bd['binnedphases'][np.argmin(orb_bd['binnedmags'])]

    # here are the points you will search.  these are already phase-ordered.
    x, y = orb_bd['binnedphases'], orb_bd['binnedmags']

    # number of points in one full cycle
    N = int(1/binsize_phase_units)

    peaks, properties = find_peaks(
        -y, prominence=prominence, width=width
    )

    # drop duplicate peaks from the wrap.  this was just to ensure you get dip
    # minima at phase=0.
    sel = (peaks <= N)

    peaks = peaks[sel]
    for k, v in properties.items():
        properties[k] = v[sel]
    params = "left_bases,right_bases,widths,left_ips,right_ips".split(",")
    for param in params:
        properties[param+"_phaseunits"] = properties[param]*binsize_phase_units

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
        'prominence': prominence,
        'width': width
    }

    return r


def prepare_cpv_light_curve(lcpath, cachedir):
    """
    Given a SPOC 2-minute light curve, remove non-zero quality flags,
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
    time = d['TIME']
    flux = d['PDCSAP_FLUX']
    qual = d['QUALITY']

    # remove non-zero quality flags
    sel = (qual == 0)

    x_obs = time[sel]
    y_obs = flux[sel]

    # normalize around 1
    y_obs /= np.nanmedian(y_obs)

    # NOTE: you could consider removing flares using a time-windowed slider
    # here.  however, for purposes of finding the periods, they are a small
    # enough fraction of the duty cycle that they can probably be ignored.

    # what is the cadence?
    cadence_sec = int(np.round(np.nanmedian(np.diff(x_obs))*24*60*60))

    starid = f'{ticid}_S{str(sector).zfill(4)}_{cadence_sec}sec'

    #
    # "light" detrending by default. (& cache it)
    #
    pklpath = os.path.join(cachedir, f"{starid}_dtr_lightcurve.pkl")
    if os.path.exists(pklpath):
        print(f"Found {pklpath}, loading and continuing.")
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
            print(f'Made {pklpath}')

    return (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend,
            cadence_sec, sector, starid)
