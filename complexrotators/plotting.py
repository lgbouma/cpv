"""
Contents:
    plot_TEMPLATE

    plot_river
    plot_phase
"""
import os, corner, pickle
from glob import glob
from datetime import datetime
import numpy as np, matplotlib.pyplot as plt, pandas as pd, pymc3 as pm
from numpy import array as nparr

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


def plot_phase(outdir):

    # get data (20sec, hardcopy for TIC 262400835)
    df = pd.read_csv('/Users/luke/Dropbox/proj/complexrotators/results/river/tic_262400835/HARDCOPY_spoc_tess_lightcurve_median_PDCSAP_FLUX_allsector_sigclipped.csv')
    time, flux = np.array(df.time), np.array(df.flux)

    t0 = np.nanmin(time) + 0.5/24
    period = 0.29822 # manually tuned...

    outpath = os.path.join(outdir, 'TIC262400835_phase.png')

    titlestr = 'TIC262400835 S32 20sec'

    plot_phased_light_curve(time, flux, t0, period, outpath, titlestr=titlestr)


def plot_multicolor_phase(outdir):
    raise NotImplementedError

    # get data (20sec, hardcopy for TIC 262400835)
    df0 = pd.read_csv('/Users/luke/Dropbox/proj/complexrotators/results/river/tic_262400835/HARDCOPY_spoc_tess_lightcurve_median_PDCSAP_FLUX_allsector_sigclipped.csv')
    time, flux = np.array(df.time), np.array(df.flux)

    t0 = np.nanmin(time) + 0.5/24
    period = 0.29822 # manually tuned...

    outpath = os.path.join(outdir, 'TIC262400835_phase.png')

    titlestr = 'TIC262400835 S32 20sec'

    plot_phased_light_curve(time, flux, t0, period, outpath, titlestr=titlestr)







def plot_phased_light_curve(
    time, flux, t0, period, outpath,
    ylimd=None, binsize_minutes=2, BINMS=2, titlestr=None
):

    x,y = time, flux-np.nanmean(flux)

    # make plot
    plt.close('all')
    set_style()

    fig, ax = plt.subplots(figsize=(4,3))

    # time units
    # x_fold = (x - t0 + 0.5 * period) % period - 0.5 * period
    # phase units
    _pd = phase_magseries(x, y, period, t0, wrap=True, sort=True)
    x_fold = _pd['phase']
    y = _pd['mags']

    if len(x_fold) > int(2e4):
        # see https://github.com/matplotlib/matplotlib/issues/5907
        mpl.rcParams['agg.path.chunksize'] = 10000

    #
    # begin the plot!
    #
    ax.scatter(x_fold, 1e2*(y), color="darkgray", label="data", marker='.',
               s=1, rasterized=True, alpha=0.3, linewidths=0)

    binsize_days = (binsize_minutes / (60*24))
    orb_bd = phase_bin_magseries(
        x_fold, y, binsize=binsize_days, minbinelems=3
    )
    ax.scatter(
        orb_bd['binnedphases'], 1e2*(orb_bd['binnedmags']), color='k',
        s=BINMS, linewidths=0,
        alpha=1, zorder=1002#, linewidths=0.2, edgecolors='white'
    )
    txt = f'$t_0$ [BTJD]: {t0-2457000:.6f}\n$P$: {period:.6f} d'
    ax.text(0.97,0.03,txt,
            transform=ax.transAxes,
            ha='right',va='bottom', color='k', fontsize='xx-small')


    if isinstance(titlestr,str):
        ax.set_title(titlestr, fontsize='small')
    ax.set_ylabel(r"Relative flux [$\times 10^{-2}$]")
    ax.set_xlabel("Phase")


    if isinstance(ylimd, dict):
        a.set_ylim(ylimd[k])
    format_ax(ax)

    fig.tight_layout()

    savefig(fig, outpath, dpi=350)
    plt.close('all')
