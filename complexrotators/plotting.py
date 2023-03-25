"""
Contents:

    | plot_TEMPLATE

    | plot_river
    | plot_phase
    | plot_phase_timegroups
    | plot_multicolor_phase
    | plot_phased_light_curve

    | plot_dipcountercheck
    | plot_cpvvetter
"""

#######################################
# ASTROBASE IMPORTS CAN BREAK LOGGING #
#######################################
from astrobase.lcmath import (
    phase_magseries, phase_bin_magseries, sigclip_magseries,
    find_lc_timegroups, phase_magseries_with_errs, time_bin_magseries
)

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
import os, pickle
from os.path import join
from glob import glob
from datetime import datetime
from copy import deepcopy
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from numpy import array as nparr
from collections import OrderedDict

import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

from astropy import units as u, constants as const
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table

import matplotlib.patheffects as pe
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import blended_transform_factory

from aesthetic.plot import savefig, format_ax, set_style

from complexrotators.paths import DATADIR, PHOTDIR
from complexrotators.lcprocessing import cpv_periodsearch

from cdips.lcproc import detrend as dtr

def plot_TEMPLATE(outdir):

    # get data

    # make plot
    plt.close('all')
    set_style()

    fig, ax = plt.subplots(figsize=(4,3))
    fig = plt.figure(figsize=(4,3))
    axd = fig.subplot_mosaic(
        """
        AB
        CD
        """#,
        #gridspec_kw={
        #    "width_ratios": [1, 1, 1, 1]
        #}
    )


    # set naming options
    s = ''

    bn = inspect.stack()[0][3].split("_")[1]
    outpath = join(outdir, f'{bn}{s}.png')
    savefig(fig, outpath, dpi=400)


def plot_river(time, flux, period, outdir, titlestr=None, cmap='Blues_r',
               cyclewindow=None, idstr=None):
    """
    Make a river plot
    """

    t0 = np.nanmin(time)
    flux -= np.nanmedian(flux) # already normalized

    cadence = np.nanmedian(np.diff(time)) # 2 minutes, in units of days
    N_obs_per_cycle = int(period / cadence)

    cycle_number = np.floor( (time-t0) / period)

    cycle_min = int(np.min(cycle_number))
    cycle_max = int(np.max(cycle_number))

    flux_arr = np.zeros(
        (N_obs_per_cycle, cycle_max-cycle_min)
    )

    for cycle_ind in range(cycle_min, cycle_max):

        begin = t0 + period*cycle_ind
        end = t0 + period*(cycle_ind+1)

        sel = (time > begin) & (time <= end)

        if len(flux[sel])/N_obs_per_cycle < 0.9:
            # for significantly empty cycles, do all nan.  "significantly
            # empty" here means any more than 5 cadences (10 minutes, out of a
            # ~1 day periods typically) off.
            flux_arr[:, cycle_ind] = 0

        elif len(flux[sel]) <= (N_obs_per_cycle-1):
            # for cycles a few cadences short, pad with a few nans at the end
            # of the array
            use_flux = np.ones(N_obs_per_cycle)*0

            # the beginning of the array is the flux
            use_flux[:flux[sel].shape[0]] = flux[sel]

            flux_arr[:, cycle_ind] = use_flux

        elif len(flux[sel]) > N_obs_per_cycle:
            use_flux = flux[sel][:N_obs_per_cycle]
            flux_arr[:, cycle_ind] = use_flux

        else:
            use_flux = flux[sel]
            flux_arr[:, cycle_ind] = use_flux

    vmin = np.nanmedian(flux)-4*np.nanstd(flux)
    vmax = np.nanmedian(flux)+4*np.nanstd(flux)

    fig, ax = plt.subplots(figsize=(6,10))
    c = ax.pcolor(np.arange(0, period, cadence),
                  list(range(cycle_min, cycle_max)),
                  flux_arr.T,
                  cmap=cmap, vmin=vmin, vmax=vmax,
                  shading='auto')

    divider0 = make_axes_locatable(ax)
    cax0 = divider0.append_axes('right', size='5%', pad=0.05)
    cb0 = fig.colorbar(c, ax=ax, cax=cax0, extend='both')
    cb0.set_label("Relative flux", rotation=270, labelpad=10)

    if isinstance(titlestr, str):
        ax.set_title(titlestr.replace("_", " "))
    ax.set_ylabel('Cycle number')
    ax.set_xlabel('Time [days]')

    if isinstance(cyclewindow, tuple):
        ax.set_ylim(cyclewindow)

    if isinstance(idstr, str):
        estr = ''
        if isinstance(cyclewindow, tuple):
            estr += (
                "_"+repr(cyclewindow).
                replace(', ','_').replace('(','').replace(')','')
            )
        outpath = join(outdir, f'{idstr}_river_{cmap}{estr}.png')
    else:
        raise NotImplementedError

    savefig(fig, outpath, writepdf=0)


def _get_cpv_lclist(lc_cadences, ticid):

    from complexrotators.getters import (
        get_2min_cadence_spoc_tess_lightcurve,
        get_20sec_cadence_spoc_tess_lightcurve
    )

    #
    # get the light curves for all desired cadences
    #
    lc_cadences = lc_cadences.split('_')
    lclist_2min, lclist_20sec = [], []
    for lctype in lc_cadences:
        print(f'Getting {lctype} cadence light curves...')
        lk_searchstr = ticid.replace('_', ' ')
        if lctype == '2min':
            lclist_2min = get_2min_cadence_spoc_tess_lightcurve(lk_searchstr)
        elif lctype == '20sec':
            lclist_20sec = get_20sec_cadence_spoc_tess_lightcurve(lk_searchstr)
        else:
            raise NotImplementedError
    print('Done getting light curves!')
    lclist = lclist_20sec + lclist_2min

    return lclist


def prepare_given_lightkurve_lc(lc, ticid, outdir):

    # metadata
    sector = lc.meta['SECTOR']

    # light curve data
    time = lc.time.value
    flux = lc.pdcsap_flux.value
    qual = lc.quality.value

    # remove non-zero quality flags
    sel = (qual == 0)

    x_obs = time[sel]
    y_obs = flux[sel]

    # normalize around 1
    y_obs /= np.nanmedian(y_obs)

    # what is the cadence?
    cadence_sec = int(np.round(np.nanmedian(np.diff(x_obs))*24*60*60))

    starid = f'{ticid}_S{str(sector).zfill(4)}_{cadence_sec}sec'

    #
    # "light" detrending by default. (& cache it)
    #
    pklpath = join(outdir, f"{starid}_dtr_lightcurve.pkl")
    if os.path.exists(pklpath):
        print(f"Found {pklpath}, loading and continuing.")
        with open(pklpath, 'rb') as f:
            lcd = pickle.load(f)
        y_flat = lcd['y_flat']
        y_trend = lcd['y_trend']
        x_trend = lcd['x_trend']
    else:
        y_flat, y_trend = dtr.detrend_flux(
            x_obs, y_obs, method='biweight', cval=2, window_length=4.0,
            break_tolerance=0.5
        )
        x_trend = deepcopy(x_obs)
        lcd = {
            'y_flat':y_flat,
            'y_trend':y_trend,
            'x_trend':x_trend
        }
        with open(pklpath, 'wb') as f:
            pickle.dump(lcd, f)
            print(f'Made {pklpath}')


    return (
        time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend,
        cadence_sec, sector, starid
    )


def plot_phase(
    outdir,
    ticid=None,
    lc_cadences='2min_20sec',
    manual_period=None,
    t0='binmin',
    ylim=None,
    binsize_minutes=10,
    xlim=[-0.6,0.6]
    ):
    """
    lc_cadences:
        string like "2min_20sec", "30min_2min_20sec", etc, for the types of
        light curves to pull for the given TIC ID.

    t0:
        - None defaults to 1618.
        - "binmin" defaults to phase-folding, and taking the arg-minimum
        - Any int or float will be passed as the manual phase.
    """

    lclist = _get_cpv_lclist(lc_cadences, ticid)

    if len(lclist) == 0:
        print(f'WRN! Did not find light curves for {ticid}. Escaping.')
        return 0

    # for each light curve (sector / cadence specific), detrend if needed, get
    # the best period, and then phase-fold.
    for lc in lclist:

        (time, flux, qual, x_obs, y_obs, y_flat,
         y_trend, x_trend, cadence_sec, sector,
         starid) = prepare_given_lightkurve_lc(lc, ticid, outdir)

        # get t0, period, lsp
        d = cpv_periodsearch(x_obs, y_flat, starid, outdir, t0=t0)

        # make the quicklook plot
        outpath = join(outdir, f'{starid}_quicklook.png')
        titlestr = starid.replace('_',' ')
        if not os.path.exists(outpath):
            plot_quicklook_cr(
                x_obs, y_obs, x_trend, y_trend, d['times'], d['fluxs'], outpath,
                titlestr
            )

        # make the phased plot
        outpath = join(
            outdir, f'{ticid}_S{str(sector).zfill(4)}_{cadence_sec}sec_phase.png'
        )
        t0 = d['t0']
        period = d['period']
        if isinstance(manual_period, float):
            period = manual_period

        plot_phased_light_curve(
            d['times'], d['fluxs'], t0, period, outpath,
            titlestr=titlestr, ylim=ylim, binsize_minutes=binsize_minutes,
            xlim=xlim
        )


def plot_phase_timegroups(
    outdir,
    ticid=None,
    lc_cadences='2min_20sec',
    manual_period=None,
    t0='binmin',
    ylim=None,
    binsize_minutes=10,
    xlim=[-0.6,0.6],
    yoffset=5
    ):
    """
    As in plot_phase
    """

    lclist = _get_cpv_lclist(lc_cadences, ticid)

    if len(lclist) == 0:
        print(f'WRN! Did not find light curves for {ticid}. Escaping.')
        return 0

    # for each light curve (sector / cadence specific), detrend if needed, get
    # the best period.
    _times, _fluxs, _t0s, _periods, _titlestrs = [],[],[],[], []

    for lc in lclist:

        (time, flux, qual, x_obs, y_obs, y_flat,
         y_trend, x_trend, cadence_sec, sector,
         starid) = prepare_given_lightkurve_lc(lc, ticid, outdir)

        # get t0, period, lsp
        d = cpv_periodsearch(x_obs, y_flat, starid, outdir, t0=t0)

        _t0 = d['t0']
        period = d['period']
        if isinstance(manual_period, float):
            period = manual_period
        titlestr = starid.replace('_',' ')

        _times.append(d['times'])
        _fluxs.append(d['fluxs'])
        _t0s.append(_t0)
        _periods.append(period)
        _titlestrs.append(titlestr)

    # merge lightcurve data, and split before making the plot.
    times = np.hstack(_times)
    fluxs = np.hstack(_fluxs)
    t0s = np.hstack(_t0s)
    periods = np.hstack(_periods)
    titlestrs = np.hstack(_titlestrs)

    from astrobase.lcmath import find_lc_timegroups
    ngroups, groups = find_lc_timegroups(times, mingap=3/24)

    # Make plots
    plt.close('all')
    set_style("science")
    fig, ax = plt.subplots(figsize=(3,7))

    plot_period = np.nanmean(_periods) + 0.002/24
    plot_t0 = t0s[0] + 0.25*plot_period
    plot_period_std = np.nanstd(_periods)

    ix = 0
    for n, g in enumerate(groups):

        cadence_day = 2/(60*24)
        N_cycles_in_group = len(times[g]) * cadence_day / plot_period

        gtime = times[g]
        gflux = fluxs[g]

        # t = t0 + P*e
        e_start = int(np.floor((gtime[0] - plot_t0)/plot_period))
        e_end = int(np.floor((gtime[-1] - plot_t0)/plot_period))
        txt = f"{e_start} - {e_end}"

        if N_cycles_in_group < 2:
            continue

        plot_phased_light_curve(
            gtime, gflux, plot_t0, plot_period, None,
            fig=fig, ax=ax,
            titlestr=None, binsize_minutes=binsize_minutes,
            xlim=xlim, yoffset=ix, showtext=txt, savethefigure=False,
            dy=yoffset
        )

        ix -= yoffset

    ax.set_title(
        f"{ticid.replace('_', ' ')} "
        f"P={plot_period*24:.3f}$\pm${plot_period_std*24:.3f}hr"
    )

    fig.text(-0.01,0.5, r"Relative flux [%]", va='center',
             rotation=90)

    if isinstance(ylim, list):
        ax.set_ylim(ylim)

    ax.set_xlabel('Phase')

    format_ax(ax)
    fig.tight_layout()

    outpath = join(outdir, f"{ticid}_P{plot_period*24:.3f}_phase_timegroups.png")
    savefig(fig, outpath, dpi=350)


def plot_quicklook_cr(x_obs, y_obs, x_trend, y_trend, x_flat, y_flat, outpath,
                      titlestr):

    from astrobase.lcmath import find_lc_timegroups
    ngroups, groups = find_lc_timegroups(x_obs, mingap=3/24)

    plt.close('all')
    set_style()
    fig, axs = plt.subplots(figsize=(12,7), nrows=2)

    for g in groups:

        ax = axs[0]
        ax.scatter(x_obs, 1e2*(y_obs-1), c="k", s=0.5, rasterized=True,
                   linewidths=0, zorder=42)
        ax.plot(x_trend, 1e2*(y_trend-1), color="C2", zorder=43, lw=0.5)
        ax.set_xticklabels([])

        # C: data - GP on same 100 day slice.  see residual transits IN THE DATA.
        ax = axs[1]
        ax.axhline(0, color="#aaaaaa", lw=0.5, zorder=-1)
        ax.scatter(x_flat, 1e2*(y_flat - 1), c="k", s=0.5, rasterized=True,
                   linewidths=0, zorder=42)

    fig.text(-0.01,0.5, r"Relative flux [%]", va='center',
             rotation=90)

    ax.set_xlabel('Time [TJD]')

    format_ax(ax)
    fig.tight_layout()

    savefig(fig, outpath, dpi=350)




def plot_phased_light_curve(
    time, flux, t0, period, outpath,
    ylim=None, xlim=[-0.6,0.6], binsize_minutes=2, BINMS=2, titlestr=None,
    showtext=True, showtitle=False, figsize=None,
    c0='darkgray', alpha0=0.3,
    c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
    fig=None, ax=None, savethefigure=True,
    findpeaks_result=None, ylabel=None, showxticklabels=True,
    yoffset=0, dy=5
    ):
    """
    Non-obvious args:
        binsize_minutes (float): binsize in units of minutes.
        BINMS (float): markersize for binned points.
        c0 (str): color of non-binned points.
        alpha0 (float): alpha for non-binned points.
        c1 (str): color of binned points.
        alpha1 (float): alpha for -binned points.
        phasewrap: [-1,1] or [0,1] in phase.
        plotnotscatter: if True, uses ax.plot to show non-binned points
        savethefigure: if False, returns fig and ax
        fig, ax: if passed, overrides default
        showtext: bool or stirng.  If string, the string is shown.
        findpeaks_result: output from ``count_phased_local_minima``
    """

    # make plot
    if savethefigure:
        plt.close('all')
        set_style("clean")

    if fig is None and ax is None:
        if figsize is None:
            fig, ax = plt.subplots(figsize=(4,3))
        else:
            fig, ax = plt.subplots(figsize=figsize)
    else:
        # if e.g., just an axis is passed, then plot on it
        pass

    x,y = time, flux-np.nanmean(flux)

    # time units
    # x_fold = (x - t0 + 0.5 * period) % period - 0.5 * period
    # phase units
    if plotnotscatter:
        _pd = phase_magseries(x, y, period, t0, wrap=phasewrap,
                              sort=False)
    else:
        _pd = phase_magseries(x, y, period, t0, wrap=phasewrap,
                              sort=True)
    x_fold = _pd['phase']
    y = _pd['mags']

    if len(x_fold) > int(2e4):
        # see https://github.com/matplotlib/matplotlib/issues/5907
        mpl.rcParams['agg.path.chunksize'] = 10000

    #
    # begin the plot!
    #
    if not plotnotscatter:
        ax.scatter(x_fold, 1e2*(y)+yoffset, color=c0, label="data", marker='.',
                   s=1, rasterized=True, alpha=alpha0, linewidths=0)
    else:
        ax.plot(x_fold, 1e2*(y)+yoffset, color=c0, label="data",
                lw=0.5, rasterized=True, alpha=alpha0)

    bs_days = (binsize_minutes / (60*24))
    orb_bd = phase_bin_magseries(x_fold, y, binsize=bs_days, minbinelems=3)
    ax.scatter(
        orb_bd['binnedphases'], 1e2*(orb_bd['binnedmags'])+yoffset, color=c1,
        s=BINMS, linewidths=0,
        alpha=alpha1, zorder=1002#, linewidths=0.2, edgecolors='white'
    )

    if isinstance(findpeaks_result, dict):

        if not xlim == [-0.6,0.6]:
            raise AssertionError("phasing off otherwise")

        tform = blended_transform_factory(ax.transData, ax.transAxes)
        peakloc = findpeaks_result['peaks_phaseunits']
        peakloc[peakloc>0.5] -= 1 # there must be a smarter way
        ax.scatter(peakloc, np.ones(len(peakloc))*0.98+yoffset, transform=tform,
                   marker='v', s=10, linewidths=0, edgecolors='none',
                   color='k', alpha=0.5, zorder=1000)

    if showtext:
        if isinstance(showtext, str):
            txt = showtext
            fontsize = 'xx-small'
            tform = blended_transform_factory(ax.transAxes, ax.transData)
            props = dict(boxstyle='square', facecolor='white', alpha=0.7, pad=0.15,
                         linewidth=0)
            sel = (orb_bd['binnedphases'] > 0.4) & (orb_bd['binnedphases'] < 0.5)
            ax.text(0.97,
                    np.nanmin(1e2*(orb_bd['binnedmags'][sel])+yoffset-dy/5), txt,
                    transform=tform, ha='right',va='top', color='k',
                    fontsize=fontsize, bbox=props)

        else:
            if isinstance(t0, float):
                #txt = f'$t_0$ [BTJD]: {t0:.6f}\n$P$: {period:.6f} d'
                txt = f'$P$: {period*24:.2f} hr' # simpler label
                fontsize='small'
            elif isinstance(t0, int):
                txt = f'$t_0$ [BTJD]: {t0:.1f}\n$P$: {period:.6f} d'
                fontsize='xx-small'
            ax.text(0.97,0.03,txt,
                    transform=ax.transAxes,
                    ha='right',va='bottom', color='k', fontsize=fontsize)

    if showtitle:
        txt = f'$t_0$ [BTJD]: {t0:.6f}. $P$: {period:.6f} d'
        ax.set_title(txt, fontsize='small')

    if isinstance(titlestr,str):
        ax.set_title(titlestr.replace("_"," "), fontsize='small')

    if savethefigure:
        ax.set_ylabel(r"Flux [$\times 10^{-2}$]")
        ax.set_xlabel("Phase")

    ax.set_xlim(xlim)
    if xlim == [-0.6,0.6]:
        ax.set_xticks([-0.5, -0.25, 0, 0.25, 0.5])
        ax.set_xticklabels([-0.5, -0.25, 0, 0.25, 0.5])

    if isinstance(ylim, (list, tuple)):
        ax.set_ylim(ylim)
    else:
        if savethefigure:
            ylim = get_ylimguess(1e2*orb_bd['binnedmags'])
            ax.set_ylim(ylim)
        else:
            pass

    if not showxticklabels:
        ax.set_xticklabels([])

    if fig is not None:
        fig.tight_layout()

    if savethefigure:
        savefig(fig, outpath, dpi=350)
        plt.close('all')
    else:
        if fig is not None:
            return fig, ax
        else:
            pass


def get_ylimguess(y):
    ylow = np.nanpercentile(y, 0.1)
    yhigh = np.nanpercentile(y, 99.9)
    ydiff = (yhigh-ylow)
    ymin = ylow - 0.35*ydiff
    ymax = yhigh + 0.35*ydiff
    return [ymin,ymax]


def _get_all_tic262400835_data():
    # get median-filtered tess 20 second data. (4.5day blended rotation
    # removed).
    df0 = pd.read_csv(
        '/Users/luke/Dropbox/proj/complexrotators/results/river/tic_262400835/'
        'HARDCOPY_spoc_tess_lightcurve_median_PDCSAP_FLUX_allsector_sigclipped.csv'
    )
    time, flux = np.array(df0.time), np.array(df0.flux)
    tess_texp = np.nanmedian(np.diff(time))

    datasets = OrderedDict()
    datasets['tess20sec'] = [time, flux, None, tess_texp]

    # get MUSCAT1 data
    m1_files = glob(
        join(PHOTDIR, 'TIC262400835_*muscat_*.csv')
    )
    for m1_file in m1_files:
        color = os.path.basename(m1_file).split("_")[3]
        timestamp = os.path.basename(m1_file).split("_")[1]
        assert color in ['g','r','z']
        datasetname = 'M1_'+color+"_"+timestamp
        _df = pd.read_csv(m1_file)
        time, flux, err = nparr(_df.BJD_TDB), nparr(_df.Flux), nparr(_df.Err)
        texp = np.nanmedian(np.diff(time))
        datasets[datasetname] = [time, flux, err, texp]

    # get MUSCAT2 data
    m2_files = glob(
        join(PHOTDIR, 'tic262400835_*achromatic*.fits')
    )
    for m2_file in m2_files:
        timestamp = os.path.basename(m2_file).split("_")[1]
        hl = fits.open(m2_file)
        # two transits...
        for color in ['g','r','z_s']:
            fitskey = 'flux_'+color
            datasetname = 'M2_'+color.replace('_s','')+"_"+timestamp
            time, flux, err = (
                nparr(hl[fitskey].data['time_bjd']),
                nparr(hl[fitskey].data['flux']),
                None
            )
            texp = np.nanmedian(np.diff(time))
            datasets[datasetname] = [time, flux, err, texp]

    return datasets


def plot_multicolor_phase(outdir, BINMS=2):
    # currently hard-coded for TIC 262400835

    #keys: tess20sec, M1_g|r|z_201216, M2_g|r|z_201213|201215
    datasets = _get_all_tic262400835_data()

    t0 = np.nanmin(datasets['tess20sec'][0]) + 0.5/24
    period = 0.29822 # manually tuned...

    # vision:
    # TESS phase-fold (binned and raw)
    # g-r phase-fold
    # g-z phase-fold
    # ... but first, just look at the data.

    # initialize
    outpath = join(outdir, 'TIC262400835_multicolor_phase_stacked.png')
    time, flux = datasets['tess20sec'][0], datasets['tess20sec'][1]
    fig, ax = plot_phased_light_curve(time, flux, t0, period, None,
                                      savethefigure=False, figsize=(4,8),
                                      showtext=False)
    keylist = [
        'M2_g_201213', 'M2_r_201213', 'M2_z_201213',
        'M2_g_201215', 'M2_r_201215', 'M2_z_201215',
        'M1_g_201216', 'M1_r_201216', 'M1_z_201216'
    ]
    y_offset = -10
    BPTOCOLORDICT = {
        'g': 'b', 'r': 'g', 'z': 'r'
    }
    ix = 0
    orb_bds = {}
    for k in keylist:
        bandpass = k.split("_")[1]
        color = BPTOCOLORDICT[bandpass]
        time, flux = datasets[k][0], datasets[k][1]
        x,y = time, flux-np.nanmean(flux)
        _pd = phase_magseries(x, y, period, t0, wrap=True, sort=True)
        x_fold, y = _pd['phase'], _pd['mags']
        ax.scatter(x_fold, 1e2*(y)+y_offset, color=color, label="data", marker='.',
                   s=1, rasterized=True, alpha=0.3, linewidths=0)

        binsize_minutes = 15
        bs_days = (binsize_minutes / (60*24))
        orb_bd = phase_bin_magseries(x_fold, y, binsize=bs_days, minbinelems=1)
        ax.scatter(
            orb_bd['binnedphases'], 1e2*(orb_bd['binnedmags'])+y_offset, color=color,
            s=BINMS, linewidths=0, alpha=1, zorder=1002
        )
        orb_bds[k] = orb_bd

        props = dict(boxstyle='square', facecolor='white', alpha=0.7, pad=0.15,
                     linewidth=0)
        txt = k
        ax.text(.95, 2.5+np.nanmedian(1e2*(y)+y_offset), txt,
                ha='right',va='center', color='k',
                fontsize='xx-small', bbox=props, zorder=1003)

        if ix % 3 != 2:
            y_offset -= 5
        else:
            y_offset -= 10

        ix += 1

    savefig(fig, outpath, dpi=350)
    plt.close('all')

    # now do it in bonafide color
    outpath = join(outdir, 'TIC262400835_multicolor_phase.png')

    time, flux = datasets['tess20sec'][0], datasets['tess20sec'][1]
    fig, ax = plot_phased_light_curve(time, flux, t0, period, None,
                                      savethefigure=False, figsize=(4,8),
                                      showtext=False)
    keylist = [
        ('M2_g_201213','M2_r_201213'),
        ('M2_r_201213','M2_z_201213'),
        ('M2_g_201213','M2_z_201213'),
        ('M2_g_201215','M2_r_201215'),
        ('M2_r_201215','M2_z_201215'),
        ('M2_g_201215','M2_z_201215'),
        ('M1_g_201216','M1_r_201216'),
        ('M1_r_201216','M1_z_201216'),
        ('M1_g_201216','M1_z_201216'),
    ]
    y_offset = -10
    ix = 0
    for b0, b1 in keylist:

        c0, c1 = b0.split("_")[1], b1.split("_")[1]
        x0, y0 = orb_bds[b0]['binnedphases'], orb_bds[b0]['binnedmags']
        x1, y1 = orb_bds[b1]['binnedphases'], orb_bds[b1]['binnedmags']

        # average of phases...
        x = x0 + (x1-x0)/2
        # the ~color
        y = ( (1+y0) / (1+y1) ) - 1

        ax.scatter(
            x, 1e2*(y)+y_offset, color='k',
            s=BINMS, linewidths=0, alpha=1, zorder=1002
        )

        props = dict(boxstyle='square', facecolor='white', alpha=0.7, pad=0.15,
                     linewidth=0)
        txt = b0.split("_")[0] + "_" + b0.split("_")[-1] + "_"+ c0 + "/" + c1
        ax.text(.95, 2.5+np.nanmedian(1e2*(y)+y_offset), txt,
                ha='right',va='center', color='k',
                fontsize='xx-small', bbox=props, zorder=1003)

        y_offset -= 6

        ix += 1

    savefig(fig, outpath, dpi=350)
    plt.close('all')


def plot_dipcountercheck(findpeaks_result, d, eval_dict, outdir, starid):
    # args: results from count_phased_local_minima and cpv_periodsearch,
    # respectively
    # eval_dict can be None

    r = findpeaks_result
    xlim = [-0.6, 0.6]

    # get data
    outpath = join(outdir, f'{starid}_dipcountercheck.png')
    if os.path.exists(outpath):
        print(f'found {outpath}, skipping')
        return

    # make plot
    plt.close('all')
    set_style('clean')

    fig = plt.figure(figsize=(4,8))
    axd = fig.subplot_mosaic(
        """
        A
        B
        C
        """
    )

    # the data
    ax = axd['A']

    ax.scatter(
        r['binned_phase'], r['binned_orig_flux'], c='k', zorder=1, s=3
    )
    ax.plot(
        r['binned_phase'], r['binned_trend_flux'], c='darkgray', zorder=2
    )
    ax.set_ylabel('Normalized Flux')

    txt = f'$P$: {d["period"]*24:.2f} hr' # simpler label
    ax.text(0.97,0.03,txt,
            transform=ax.transAxes,
            ha='right',va='bottom', color='k', fontsize='small')

    ax.set_title(starid.replace("_"," "), fontsize='small')

    # the transformation, used to count dips
    ax = axd['B']
    ax.scatter(
        r['binned_phase'], r['binned_search_flux'], c='k', s=3
    )
    ax.set_ylabel("Residual")

    # show the dips
    for ax in [axd['A'], axd['B']]:
        if not xlim == [-0.6,0.6]:
            raise AssertionError("phasing off otherwise")
        ax.set_xlim(xlim)
        tform = blended_transform_factory(ax.transData, ax.transAxes)
        peakloc = r['peaks_phaseunits']
        peakloc[peakloc>0.5] -= 1 # there must be a smarter way
        ax.scatter(peakloc, np.ones(len(peakloc))*0.98, transform=tform,
                   marker='v', s=10, linewidths=0, edgecolors='none',
                   color='k', alpha=0.5, zorder=1000)

    # text summary
    ax = axd['C']
    ax.set_axis_off()
    txt0 = (
        f'{starid}\n'
        f'N dips: {r["N_peaks"]}\n'
        f'p2p: {r["p2p_est"]:.2e}, height: {r["height"]:.2e}\n'
        f'locs: {r["peaks_phaseunits"]}\n'
        f'left bases: {r["properties"]["left_bases_phaseunits"]}\n'
        f'right bases: {r["properties"]["right_bases_phaseunits"]}'
    )
    ax.text(0.03, 0.03, txt0,
            transform=ax.transAxes, ha='left', va='bottom', color='k')

    if isinstance(eval_dict, dict):
        txt1 = ''
        if eval_dict['found_correct_ndips']:
            txt1 += 'found correct ndips'
            color = 'green'
        else:
            txt1 += 'found incorrect ndips'
            color = 'red'
        ax.text(0.03, 0.97, txt1,
                transform=ax.transAxes, ha='left', va='top', color=color)

        txt2 = ''
        if eval_dict['found_toofew_dips'] or eval_dict['found_toomany_dips']:
            if eval_dict['found_toofew_dips']:
                txt2 += 'found too few dips'
            if eval_dict['found_toomany_dips']:
                txt2 += 'found too many dips'
            color = 'red'
            ax.text(0.03, 0.8, txt2,
                    transform=ax.transAxes, ha='left', va='top', color=color)

    # finish
    savefig(fig, outpath, dpi=400)


def plot_cpvvetter(
    outpath,
    lcpath,
    starid,
    periodsearch_result=None,
    findpeaks_result=None,
    binsize_minutes=10
    ):

    # get data
    d = periodsearch_result
    r = findpeaks_result
    hdul = fits.open(lcpath)
    data = hdul[1].data
    hdul.close()
    qual = data["QUALITY"]
    sel = (qual == 0)
    import IPython; IPython.embed()
    #FIXME FIXME FIXME TODO TODO TODO PULL THE DATA NEEDED TO MAKE REMAINING
    #PLOTS
    #xc, yc = 

    # make plot
    plt.close('all')
    set_style("clean")

    fig = plt.figure(figsize=(8,4.5))
    axd = fig.subplot_mosaic(
        """
        AAAABBCC
        AAAABBCC
        AAAADDDD
        EEFFKKLL
        GGHHKKLL
        IIJJKKLL
        """
    )

    # pdcsap flux vs time (qual==0, after 5-day median smooth)
    ax = axd['D']
    bd = time_bin_magseries(d['times'], d['fluxs'], binsize=600, minbinelems=1)
    ax.scatter(d['times'], d['fluxs'], c='lightgray', s=1, zorder=1)
    ax.scatter(bd['binnedtimes'], bd['binnedmags'], c='k', s=1, zorder=2)
    txt = 'PDCSAP, 5d median filter'
    bbox = dict(facecolor='white', alpha=0.9, pad=0, edgecolor='white')
    ax.text(0.03, 0.97, txt, ha='left', va='top', bbox=bbox, zorder=3,
            transform=ax.transAxes)
    ax.update({'xlabel': 'Time [BTJD]', 'ylabel': 'Flux'})

    # pdm periodogram
    ax = axd['B']
    ax.plot(d['lsp']['periods'], d['lsp']['lspvals'], c='k', lw=1)
    ax.scatter(d['lsp']['nbestperiods'][:5], d['lsp']['nbestlspvals'][:5],
               marker='v', s=5, linewidths=0, edgecolors='none',
               color='k', alpha=0.5, zorder=1000)
    ymin, ymax = ax.get_ylim()
    ax.vlines(d['period'], ymin, ymax, colors='darkgray', alpha=1,
              linestyles='-', zorder=-2, linewidths=1)
    ax.vlines([0.5*d['period'], 2*d['period']], ymin, ymax, colors='darkgray',
              alpha=0.5, linestyles=':', zorder=-2, linewidths=1)
    ax.set_ylim([ymin, ymax])
    ax.update({'xlabel': 'Period [d]', 'ylabel': 'PDM Θ', 'xscale': 'log'})

    # phased LC
    ax = axd['A']
    plot_phased_light_curve(
        d['times'], d['fluxs'], d['t0'], d['period'], None, ylim=None,
        xlim=[-0.6,0.6], binsize_minutes=2, BINMS=2, titlestr=None,
        showtext=True, showtitle=False, figsize=None, c0='darkgray',
        alpha0=0.3, c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
        fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
        showxticklabels=False
    )

    # phased LC at 2x period
    ax = axd['F']
    plot_phased_light_curve(
        d['times'], d['fluxs'], d['t0'], 2*d['period'], None, ylim=None,
        xlim=[-0.6,0.6], binsize_minutes=2, BINMS=2, titlestr=None,
        showtext=True, showtitle=False, figsize=None, c0='darkgray',
        alpha0=0.3, c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
        fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
        showxticklabels=False, ylabel='Flux'
    )
    ax.set_xlabel("")

    # phased LC at 0.5x period
    ax = axd['H']
    plot_phased_light_curve(
        d['times'], d['fluxs'], d['t0'], 0.5*d['period'], None, ylim=None,
        xlim=[-0.6,0.6], binsize_minutes=2, BINMS=2, titlestr=None,
        showtext=True, showtitle=False, figsize=None, c0='darkgray',
        alpha0=0.3, c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
        fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
        showxticklabels=False, ylabel='Flux'
    )
    ax.set_xlabel("")

    # # phased norm LC, w/ smooth model
    # ax = axd['E']

    # # resids of above, w/ ticks for auto dip find
    # ax = axd['G']

    # # phased flux in BGD aperture
    ax = axd['I']
    plot_phased_light_curve(
        d['times'], d['fluxs'], d['t0'], d['period'], None, ylim=None,
        xlim=[-0.6,0.6], binsize_minutes=2, BINMS=2, titlestr=None,
        showtext=True, showtitle=False, figsize=None, c0='darkgray',
        alpha0=0.3, c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
        fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
        showxticklabels=False
    )

    # # phased ctd x/y
    # ax = axd['J']

    # # star info (Gaia, TIC8, dip search)
    # ax = axd['L']

    # # gaia CAMD
    # ax = axd['K']

    # # DSS query
    # ax = axd['C']


    # set naming options
    s = ''

    fig.tight_layout()

    savefig(fig, outpath)
