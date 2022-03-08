"""
Contents:
    plot_TEMPLATE

    plot_river
    plot_phase
    plot_multicolor_phase
    plot_phased_light_curve
"""
import os, pickle
from glob import glob
from datetime import datetime
from copy import deepcopy
import numpy as np, matplotlib.pyplot as plt, pandas as pd, pymc3 as pm
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

from aesthetic.plot import savefig, format_ax, set_style
from astrobase.lcmath import (
    phase_magseries, phase_bin_magseries, sigclip_magseries,
    find_lc_timegroups, phase_magseries_with_errs, time_bin_magseries
)

from complexrotators.paths import DATADIR, PHOTDIR
from complexrotators.helpers import (
    get_complexrot_data,
    get_complexrot_twentysec_data
)

from complexrotators.getters import (
    get_2min_cadence_spoc_tess_lightcurve,
    get_20sec_cadence_spoc_tess_lightcurve
)
from complexrotators.lcprocessing import (
    cr_periodsearch
)

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
    outpath = os.path.join(outdir, f'{bn}{s}.png')
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

        if len(flux[sel]) < (N_obs_per_cycle-5):
            # for significantly empty cycles, do all nan.  "significantly
            # empty" here means any more than 5 cadences (10 minutes, out of a
            # ~1 day periods typically) off.
            flux_arr[:, cycle_ind] = 0

        elif len(flux[sel]) < N_obs_per_cycle:
            # for cycles a few cadences short, pad with a few nans at the end
            # of the array
            use_flux = np.ones(N_obs_per_cycle)*0

            # the beginning of the array is the flux
            use_flux[:flux[sel].shape[0]] = flux[sel]

            flux_arr[:, cycle_ind] = use_flux

        elif len(flux[sel]) < (N_obs_per_cycle+5):
            use_flux = flux[sel][:N_obs_per_cycle]
            flux_arr[:, cycle_ind] = use_flux

        elif len(flux[sel]) > (N_obs_per_cycle+5):
            raise NotImplementedError('How did this happen?')

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

    if isinstance(titlestr, str):
        ax.set_title(titlestr)
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
        outpath = os.path.join(outdir, f'{idstr}_river_{cmap}{estr}.png')
    else:
        raise NotImplementedError

    savefig(fig, outpath, writepdf=0)


def plot_phase(
    outdir,
    ticid=None,
    lc_cadences='2min_20sec',
    MANUAL_PERIOD=None,
    ylim=None,
    binsize_minutes=10
):
    """
    lc_cadences: string like "2min_20sec", "30min_2min_20sec", etc, for the types
    of light curves to pull for the given TIC ID.
    """

    # get the light curves for all desired cadences
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

    if len(lclist) == 0:
        print(f'WRN! Did not find light curves for {ticid}. Escaping.')
        return 0

    # for each light curve (sector / cadence specific), detrend if needed, get
    # the best period, and then phase-fold.
    for lc in lclist:

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

        # "light" detrending by default
        y_flat, y_trend = dtr.detrend_flux(
            x_obs, y_obs, method='biweight', cval=2, window_length=4.0,
            break_tolerance=0.5
        )
        x_trend = deepcopy(x_obs)

        # get t0, period, lsp
        d = cr_periodsearch(
            x_obs, y_flat, f'{ticid}_S{sector}_{cadence_sec}sec', outdir
        )

        starid = f'{ticid}'

        # make the quicklook plot
        outpath = os.path.join(
            outdir, f'{starid}_S{str(sector).zfill(4)}_{cadence_sec}sec_quicklook.png'
        )
        titlestr = f'{starid} S{str(sector).zfill(4)} {cadence_sec}sec'

        plot_quicklook_cr(
            x_obs, y_obs, x_trend, y_trend, d['times'], d['fluxs'], outpath,
            titlestr
        )

        # make the phased plot
        outpath = os.path.join(
            outdir, f'{starid}_S{str(sector).zfill(4)}_{cadence_sec}sec_phase.png'
        )
        period = d['period']
        if isinstance(MANUAL_PERIOD, float):
            period = MANUAL_PERIOD

        plot_phased_light_curve(
            d['times'], d['fluxs'], d['t0'], period, outpath,
            titlestr=titlestr, ylim=ylim, binsize_minutes=binsize_minutes
        )


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

    fig.text(-0.01,0.5, r"Relative flux [$\times 10^{-2}$]", va='center',
             rotation=90)

    ax.set_xlabel('Time [TJD]')

    format_ax(ax)
    fig.tight_layout()

    savefig(fig, outpath, dpi=350)




def plot_phased_light_curve(
    time, flux, t0, period, outpath,
    ylim=None, binsize_minutes=2, BINMS=2, titlestr=None,
    showtext=True, showtitle=False, figsize=None,
    c0='darkgray', alpha0=0.3,
    c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
    fig=None, ax=None, savethefigure=True
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
    """

    # make plot
    plt.close('all')
    set_style()

    if fig is None and ax is None:
        if figsize is None:
            fig, ax = plt.subplots(figsize=(4,3))
        else:
            fig, ax = plt.subplots(figsize=figsize)
    else:
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
        ax.scatter(x_fold, 1e2*(y), color=c0, label="data", marker='.',
                   s=1, rasterized=True, alpha=alpha0, linewidths=0)
    else:
        ax.plot(x_fold, 1e2*(y), color=c0, label="data",
                lw=0.5, rasterized=True, alpha=alpha0)

    bs_days = (binsize_minutes / (60*24))
    orb_bd = phase_bin_magseries(x_fold, y, binsize=bs_days, minbinelems=3)
    ax.scatter(
        orb_bd['binnedphases'], 1e2*(orb_bd['binnedmags']), color=c1,
        s=BINMS, linewidths=0,
        alpha=alpha1, zorder=1002#, linewidths=0.2, edgecolors='white'
    )
    if showtext:
        txt = f'$t_0$ [BTJD]: {t0-2457000:.6f}\n$P$: {period:.6f} d'
        ax.text(0.97,0.03,txt,
                transform=ax.transAxes,
                ha='right',va='bottom', color='k', fontsize='xx-small')
    if showtitle:
        txt = f'$t_0$ [BTJD]: {t0-2457000:.6f}. $P$: {period:.6f} d'
        ax.set_title(txt, fontsize='small')


    if isinstance(titlestr,str):
        ax.set_title(titlestr, fontsize='small')
    ax.set_ylabel(r"Relative flux [$\times 10^{-2}$]")
    ax.set_xlabel("Phase")


    if isinstance(ylim, (list, tuple)):
        ax.set_ylim(ylim)
    else:
        ylim = get_ylimguess(1e2*orb_bd['binnedmags'])
        ax.set_ylim(ylim)

    format_ax(ax)

    fig.tight_layout()

    if savethefigure:
        savefig(fig, outpath, dpi=350)
        plt.close('all')
    else:
        return fig, ax


def get_ylimguess(y):
    ylow = np.nanpercentile(y, 0.1)
    yhigh = np.nanpercentile(y, 99.9)
    ydiff = (yhigh-ylow)
    ymin = ylow - 0.3*ydiff
    ymax = yhigh + 0.3*ydiff
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
        os.path.join(PHOTDIR, 'TIC262400835_*muscat_*.csv')
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
        os.path.join(PHOTDIR, 'tic262400835_*achromatic*.fits')
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
    outpath = os.path.join(outdir, 'TIC262400835_multicolor_phase_stacked.png')
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
    outpath = os.path.join(outdir, 'TIC262400835_multicolor_phase.png')

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
