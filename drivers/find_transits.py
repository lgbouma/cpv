"""
Once you have a CPV, fit out the "mean signal", subtract, and look for transits
in the residual.
"""

#############
## LOGGING ##
#############
import logging
from complexrotators import log_sub, log_fmt, log_date_fmt

LOCAL_DEBUG = 1
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
import os, pickle
from os.path import join
from glob import glob
import numpy as np, pandas as pd, matplotlib.pyplot as plt

from complexrotators.paths import LOCALDIR, SPOCDIR, TABLEDIR
from complexrotators.getters import (
    _get_lcpaths_given_ticid, _get_local_lcpaths_given_ticid,
    _get_lcpaths_fromlightkurve_given_ticid
)
from complexrotators.lcprocessing import (
    cpv_periodsearch, count_phased_local_minima, prepare_cpv_light_curve
)
from complexrotators.plotting import (
    plot_quasiperiodic_removal_diagnostic
)

from complexrotators import pipeline_utils as pu


def find_transits(ticid, sample_id):

    cachedir = join(LOCALDIR, "cpv_transit_finding")
    if not os.path.exists(cachedir): os.mkdir(cachedir)

    cachedir = join(cachedir, sample_id)
    if not os.path.exists(cachedir): os.mkdir(cachedir)

    minexitcode = -1
    cand_logpaths = glob(join(cachedir, f"tess*00{ticid}-*transitfindingstatus.log"))
    foundexitcodes = []
    if len(cand_logpaths) > 0:
        for cand_logpath in cand_logpaths:
            st = pu.load_status(cand_logpath)
            if 'exitcode' in st:
                exitcode = st['exitcode']['exitcode']
                foundexitcodes.append(int(exitcode))
        if len(foundexitcodes) > 0:
            minexitcode = np.nanmin(foundexitcodes)

    MINIMUM_EXITCODE = 99
    # 1 if any kind of exit means do not rerun
    # 2 if only a periodogram or not enoigh dip exit means dont rerun
    if minexitcode >= MINIMUM_EXITCODE:
        LOGINFO(f"TIC{ticid}: found log for {ticid} with exitcode {minexitcode}. skip.")
        return 1

    #
    # get the light curves for all desired sectors and cadences
    #
    if LOCAL_DEBUG:
        lcpaths = _get_lcpaths_fromlightkurve_given_ticid(ticid)
    else:
        lcpaths = _get_local_lcpaths_given_ticid(ticid)

    #
    # for each light curve (sector / cadence specific), detrend, ((remove
    # flares)), get the best period, and then phase-fold.
    #
    for lcpath in lcpaths:

        # instantiate the log
        lcpbase = os.path.basename(lcpath).replace(".fits", "")
        logpath = join(cachedir, f'{lcpbase}_runstatus.log')
        if not os.path.exists(logpath):
            lcpd = {
                'lcpath': lcpath,
                'ticid': ticid
            }
            pu.save_status(logpath, 'lcpath', lcpd)
            LOGINFO(f"Made {logpath}")

        st = pu.load_status(logpath)
        if 'exitcode' in st:
            exitcode = st['exitcode']['exitcode']
            if minexitcode >= MINIMUM_EXITCODE:
                LOGINFO(f"{lcpbase}: found exitcode {exitcode}. skip.")
                continue

        # get the relevant light curve data
        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid) = prepare_cpv_light_curve(lcpath, cachedir)

        # get period, t0, and periodogram (PDM or LombScargle)
        d = cpv_periodsearch(
            x_obs, y_flat, starid, cachedir, t0='binmin', periodogram_method='pdm'
        )

        method = 'mean'
        cachepath = join(cachedir, f'{starid}_quasiperiodicremoval_{method}.pkl')
        r = remove_quasiperiodic_signal(
            d['times'], d['fluxs'], d['t0'], d['period'], cachepath,
            starid, method=method, pgdict=d
        )

        # TODO FIXME: now find planets in the residual


def remove_quasiperiodic_signal(
    time, flux, t0, period, cachepath, starid,
    method=None, make_diagnostic_plot=1, pgdict=None
):
    """
    Given np.ndarray time/flux, the period and phase, and a `cachepath` to
    store the output

    Args:

        time/flux (np.ndarray): time and flux

        t0/period (float): reference epoch and period (days)

        cachepath (str): path where a pickle file saving the contents will be
        cached.

        method (str): "mean", or "locor".  The former gets the average signal
        shape in 5-minute bins... and subtracts that.  The latter runs LOCOR.

        make_diagnostic_plot (bool): whether to make a plot that tells you how
        you did.

    Returns:
        dictionary containing keys to the subtracted (actually, divided out)
        data.
    """

    if method == 'locor':
        raise NotImplementedError

    from astrobase.lcmath import (
        phase_magseries, phase_bin_magseries, sigclip_magseries,
        find_lc_timegroups, phase_magseries_with_errs, time_bin_magseries
    )

    sel = np.isfinite(time) & np.isfinite(flux)
    time = time[sel]
    flux = flux[sel]

    x = time
    y = flux

    _pd_sw = phase_magseries(x, y, period, t0, wrap=1, sort=True)
    x_sw = _pd_sw['phase']
    y_sw = _pd_sw['mags']
    _pd_nsnw = phase_magseries(x, y, period, t0, wrap=0, sort=False)
    x_nsnw = _pd_nsnw['phase']
    y_nsnw = _pd_nsnw['mags']

    cyclenum_nsnw = np.floor((x - t0) / period)

    sortorder = np.argsort(_pd_nsnw['phase'])
    time_sw = np.concatenate((x[sortorder], x[sortorder]))
    cyclenum_sw = np.floor( (time_sw - t0) / period)

    binsize_minutes = 10
    bs_days = (binsize_minutes / (60*24))
    bd = phase_bin_magseries(x_sw, y_sw, binsize=bs_days, minbinelems=3)

    x_sw_b = bd['binnedphases']
    y_sw_b = bd['binnedmags']

    # now: get a model and interpolate!
    from scipy.interpolate import interp1d, PchipInterpolator
    interp_methods = ["linear", "slinear", "quadratic", "pchip"]

    fndict = {}
    resids = []
    for interp_method in interp_methods:
        if interp_method in ["linear", "slinear", "quadratic"]:
            fn = interp1d(x_sw_b, y_sw_b, kind=interp_method,
                          fill_value="extrapolate")
        elif interp_method == "pchip":
            fn = PchipInterpolator(x_sw_b, y_sw_b)

        fndict[interp_method] = fn

        this_y_sw_model = fn(x_sw)
        resids.append( np.sum(np.abs(this_y_sw_model - y_sw)**2) )

    best_interp_ix = np.argmin(np.array(resids))
    best_interp_key = interp_methods[best_interp_ix]

    # adopted interpolation model
    fn = fndict[best_interp_key]

    # *_nsnw: non-sorted, no wrap
    y_model_nsnw = fn(x_nsnw)
    y_model_sw = fn(x_sw)

    y_resid_nsnw = y_nsnw - y_model_nsnw
    y_resid_sw = y_sw - y_model_sw

    resid_bd_nsnw = phase_bin_magseries(x_nsnw, y_resid_nsnw, binsize=bs_days, minbinelems=3)
    resid_bd_sw = phase_bin_magseries(x_sw, y_resid_sw, binsize=bs_days, minbinelems=3)

    # estimate noise in the residual
    tb_nsnw = time_bin_magseries(time, y_resid_nsnw, binsize=3600)
    def mad(x):
        return np.nanmedian( np.abs(x - np.nanmean(x)) )
    mad_resid_1hr = mad(tb_nsnw['binnedmags'])

    x_resid_sw_b = resid_bd_sw['binnedphases']
    y_resid_sw_b = resid_bd_sw['binnedmags']

    x_resid_nsnw_b = resid_bd_nsnw['binnedphases']
    y_resid_nsnw_b = resid_bd_nsnw['binnedmags']

    x_w_model = np.linspace(-1,1,1000)
    y_w_model = fn(x_w_model)

    x_nw_model = np.linspace(0,1,1000)
    y_nw_model = fn(x_nw_model)

    assert len(time) == len(y_resid_nsnw)
    cachedir = os.path.dirname(cachepath)
    resid_pgdict = cpv_periodsearch(
        time, y_resid_nsnw, starid+"_resid0", cachedir, t0='binmin',
        periodogram_method='pdm'
    )

    out_dict = {
        'pgdict': pgdict, # nominal periodogram dictionary
        'resid_pgdict': resid_pgdict, # periodogram on the RESIDUAL dictionary
        'period': period,
        'best_interp_key': best_interp_key,
        'mad_resid_1hr': mad_resid_1hr,
        'fn': fn,
        'time': time, # input time
        'flux': flux, # input flux
        'x_sw': x_sw, # sorted & wrapped phase
        'y_sw': y_sw, # sorted & wrapped flux
        'x_sw_b': x_sw_b, # "", binned
        'y_sw_b': y_sw_b, #  "", binned
        'x_nsnw': x_nsnw, # not sorted & not wrapped phase
        'y_nsnw': y_nsnw, # not sorted & not wrapped flux
        'y_model_sw': y_model_sw, # sorted & wrapped model flux
        'y_model_nsnw': y_model_nsnw, # not sorted & not wrapped model flux
        'y_resid_sw': y_resid_sw, # y_sw - y_model_sw
        'y_resid_nsnw': y_resid_nsnw, # y_nsnw - y_model_nsnw
        'x_resid_sw_b': x_resid_sw_b,
        'y_resid_sw_b': y_resid_sw_b,
        'x_resid_nsnw_b': x_resid_nsnw_b,
        'y_resid_nsnw_b': y_resid_nsnw_b,
        'x_w_model': x_w_model,
        'y_w_model': y_w_model,
        'x_nw_model': x_nw_model,
        'y_nw_model': y_nw_model,
        'cyclenum_nsnw': cyclenum_nsnw,
        'cyclenum_sw': cyclenum_sw,
    }

    if not os.path.exists(cachepath):
        with open(cachepath, "wb") as f:
            pickle.dump(out_dict, f)
            LOGINFO(f"Wrote {cachepath}")
    else:
        LOGINFO(f"Found {cachepath}")

    if make_diagnostic_plot:
        pngpath = cachepath.replace(".pkl", "_diagnostic.png")
        plot_quasiperiodic_removal_diagnostic(out_dict, pngpath)

    return out_dict


def main():

    #sample_id = 'debug'
    #ticids = ["402980664", "224283342", "254612758"]

    sample_id = '20230613_LGB_RJ_selfnapplied'
    csvpath = join(
        TABLEDIR, "2023_catalog_table",
        "20230613_LGB_RJ_CPV_TABLE_supplemental_selfnapplied_BACKUP.csv"
    )
    df = pd.read_csv(csvpath, sep="|")
    ticids = np.array(df.ticid).astype(str)

    sample_id = 'debug'
    ticids = ['245902096']

    for ticid in ticids:
        LOGINFO(42*'-')
        LOGINFO(f"Beginning find_transits for {ticid}...")
        find_transits(ticid, sample_id)

if __name__ == "__main__":
    main()
