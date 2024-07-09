"""
Contents:

    Templates:
        | plot_TEMPLATE

    Light curves:
        | plot_river
        | plot_phase
        | plot_phase_timegroups
        | plot_multicolor_phase
        | plot_phased_light_curve
        | plot_lc_mosaic
        | plot_full_lcmosaic
        | plot_beforeafter_mosaic
        | plot_phase_timegroups_mosaic
        | plot_cadence_comparison

    Movies:
        | plot_movie_phase_timegroups
        | plot_movie_specriver

    Single-object:
        | plot_tic4029_segments

    Pipeline:
        | plot_dipcountercheck
        | plot_cpvvetter
        | plot_quasiperiodic_removal_diagnostic

    Spectra:
        | plot_spectrum_windows
        | plot_winered_windows
        | plot_lineevolnpanel

    Population:
        | plot_catalogscatter
        | plot_magnetic_bstar_comparison
"""

#######################################
# ASTROBASE IMPORTS CAN BREAK LOGGING #
#######################################
from astrobase.lcmath import (
    phase_magseries, phase_bin_magseries, sigclip_magseries,
    find_lc_timegroups, phase_magseries_with_errs, time_bin_magseries
)
from astrobase.services.identifiers import tic_to_gaiadr2
from astrobase.lcmath import sigclip_magseries

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
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams

from astropy import units as u, constants as const
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table

import matplotlib.patheffects as pe
from matplotlib.ticker import MaxNLocator, FixedLocator, FuncFormatter
from matplotlib.transforms import blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from matplotlib.ticker import (
    MultipleLocator, FormatStrFormatter, AutoMinorLocator
)


from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d, PchipInterpolator

from aesthetic.plot import savefig, format_ax, set_style

from complexrotators.paths import (
    DATADIR, PHOTDIR, SPOCDIR, LOCALDIR, TABLEDIR
)
from complexrotators.lcprocessing import cpv_periodsearch

from astroquery.mast import Catalogs
from cdips.lcproc import detrend as dtr
from cdips.utils.gaiaqueries import (
    given_dr2_sourceids_get_edr3_xmatch, given_source_ids_get_gaia_data
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
    outpath = join(outdir, f'{bn}{s}.png')
    savefig(fig, outpath, dpi=400)


def plot_river(time, flux, period, outdir, titlestr=None, cmap='Blues_r',
               cyclewindow=None, idstr=None, t0=None, savstr=None, vmin=None,
               vmax=None):
    """
    Make a river plot
    """

    finite_time = ~pd.isnull(time)

    time = time[finite_time]
    flux = flux[finite_time]

    if t0 is None:
        t0 = np.nanmin(time)

    flux -= np.nanmedian(flux) # already normalized

    cadence = np.nanmedian(np.diff(time)) # 2 minutes, in units of days
    N_obs_per_cycle = int(np.round(period / cadence))

    cycle_number = np.floor( (time-t0) / period)

    cycle_min = int(np.nanmin(cycle_number))
    cycle_max = int(np.nanmax(cycle_number))

    flux_arr = np.zeros(
        (N_obs_per_cycle, cycle_max-cycle_min)
    )

    for ix, cycle_ind in enumerate(range(cycle_min, cycle_max)):

        begin = t0 + period*cycle_ind
        end = t0 + period*(cycle_ind+1)

        sel = (time > begin) & (time <= end)

        ALGORITHM1 = 0 # janky
        ALGORITHM2 = 1 # interpolates to nearest flux value any desired grid time

        assert ALGORITHM1 + ALGORITHM2 == 1

        if ALGORITHM1:
            raise NotImplementedError("indexing by cycle_ind wrong")
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

        if ALGORITHM2:

            if len(time[sel])/N_obs_per_cycle < 0.2:
                # for significantly empty cycles, do all nan.  "significantly
                # empty" here means any more than 5 cadences (10 minutes, out of a
                # ~1 day periods typically) off.
                flux_arr[:, ix] = np.nan

            else:
                # select points slightly inside and outside of this cycle
                out_sel = (
                    (time > begin-0.03*period) & (time <= end+0.03*period)
                )

                t_start_this_cycle = t0 + cycle_ind * period
                t_end_this_cycle = t0 + (cycle_ind+1) * period

                fn = interp1d(time[out_sel], flux[out_sel],
                              kind='nearest', fill_value=np.nan,
                              bounds_error=False)

                t_grid_this_cycle = np.arange(
                    t_start_this_cycle, t_end_this_cycle, cadence
                )

                assert len(t_grid_this_cycle) == N_obs_per_cycle

                flux_this_cycle = fn(t_grid_this_cycle)

                print(cycle_ind, N_obs_per_cycle)
                N_uniq = len(np.unique(flux_this_cycle))
                N_pt = len(flux_this_cycle)

                if N_uniq / N_pt < 0.8:
                    flux_arr[:, ix] = np.nan
                else:
                    try:
                        flux_arr[:, ix] = flux_this_cycle
                    except IndexError:
                        import IPython; IPython.embed()

    if not isinstance(vmin, float):
        vmin = np.nanmedian(flux)-4*np.nanstd(flux)
    if not isinstance(vmax, float):
        vmax = np.nanmedian(flux)+4*np.nanstd(flux)

    set_style('science')
    fig, ax = plt.subplots(figsize=(6,5))
    c = ax.pcolor(np.arange(0, period, cadence)/period - 0.5,
                  list(range(cycle_min, cycle_max)),
                  100*flux_arr.T,
                  cmap=cmap,
                  vmin=vmin, vmax=vmax,
                  #norm=colors.SymLogNorm(  # NOTE: tried this; looks bad
                  #    linthresh=0.003, linscale=1,
                  #     vmin=vmin, vmax=vmax, base=10
                  #),
                  shading='auto', rasterized=True)

    divider0 = make_axes_locatable(ax)
    if isinstance(titlestr, str):
        cax0 = divider0.append_axes('right', size='5%', pad=0.05)
        cb0 = fig.colorbar(c, ax=ax, cax=cax0, extend='both')
        cb0.set_label("$\Delta$ Flux [%]", rotation=270, labelpad=10)
    elif titlestr is None:
        # sick inset colorbar
        x0,y0,dx,dy = 0.02, -0.09, 0.3, 0.02
        axins1 = inset_axes(ax, width="100%", height="100%",
                            bbox_to_anchor=(x0,y0,dx,dy),
                            loc='lower left',
                            bbox_transform=ax.transAxes)

        cb = fig.colorbar(c, cax=axins1, orientation="horizontal",
                          extend="both")
        cb.set_label("$\Delta$ Flux [%]", rotation=0, labelpad=3)
        #cb.set_ticks([min_cycle, max_cycle])
        #cb.set_ticklabels([int(min_cycle), int(max_cycle)])
        #cb.ax.tick_params(labelsize='x-small')
        #cb.ax.tick_params(size=0, which='both') # remove the ticks
        #cb.ax.yaxis.set_ticks_position('left')
        #cb.ax.yaxis.set_label_position('left')



    if isinstance(titlestr, str):
        ax.set_title(titlestr.replace("_", " "))
    ax.set_ylabel('Cycle number')
    ax.set_xlabel('Phase, $\phi$')

    if isinstance(cyclewindow, tuple):
        ax.set_ylim(cyclewindow)

    s = ''
    if isinstance(savstr, str):
        s += f"_{savstr}"
    if isinstance(vmin, float):
        s += f"_vmin{vmin:.5f}"
    if isinstance(vmax, float):
        s += f"_vmax{vmax:.5f}"
    if isinstance(cmap, str):
        cmapstr = cmap
    else:
        cmapstr = 'gist_stern_truncated'

    if isinstance(idstr, str):
        estr = ''
        if isinstance(cyclewindow, tuple):
            estr += (
                "_"+repr(cyclewindow).
                replace(', ','_').replace('(','').replace(')','')
            )
        outpath = join(outdir, f'{idstr}_river_{cmapstr}{estr}{s}.png')
    else:
        raise NotImplementedError

    fig.tight_layout()
    savefig(fig, outpath, writepdf=1)


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

    # author
    author = lc.author

    # metadata
    sector = lc.meta['SECTOR']

    # light curve data
    time = lc.time.value
    if "SPOC" in lc.author:
        flux = lc.pdcsap_flux.value
    elif lc.author == "QLP":
        flux = lc.flux.value
    else:
        raise NotImplementedError
    qual = lc.quality.value

    # remove non-zero quality flags
    sel = (qual == 0)

    x_obs = time[sel]
    y_obs = flux[sel]

    # normalize around 1
    y_obs /= np.nanmedian(y_obs)

    # what is the cadence?
    cadence_sec = int(np.round(np.nanmedian(np.diff(x_obs))*24*60*60))

    starid = f'{ticid}_S{str(sector).zfill(4)}_{author}_{cadence_sec}sec'

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
    binsize_phase=0.005,
    xlim=[-0.6,0.6],
    t0_offset=None
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

        if isinstance(t0_offset, (int,float)):
            t0 += t0_offset

        print(42*'-')
        print(t0, period)

        plt.close("all")
        plot_phased_light_curve(
            d['times'], d['fluxs'], t0, period, outpath,
            titlestr=titlestr, ylim=ylim, binsize_phase=binsize_phase,
            xlim=xlim, savethefigure=1
        )


def plot_phase_timegroups(
    outdir,
    ticid=None,
    lc_cadences='2min_20sec',
    manual_period=None,
    t0='binmin',
    ylim=None,
    binsize_phase=0.005,
    xlim=[-0.6,0.6],
    yoffset=5,
    showtitle=1,
    figsize_y=7,
    do4029_resid=0,
    cyclewindow=None,
    mingap=3/24,  # days
    rasterized=False
    ):
    """
    As in plot_phase
    """

    lclist = _get_cpv_lclist(lc_cadences, ticid)

    ## NOTE: example manual download of a particular subset
    #import lightkurve as lk
    #r = lk.search_lightcurve("TIC 300651846")
    #lcc = r[r.author == 'QLP'].download_all()
    #lclist = [_l for _l in lcc]

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
        if not isinstance(t0, float) and isinstance(manual_period, float):
            d = cpv_periodsearch(x_obs, y_flat, starid, outdir, t0=t0)
        else:
            d = {'times': x_obs, 'fluxs': y_flat,
                 't0': t0, 'period': manual_period}

        _t0 = d['t0']
        if isinstance(t0, float):
            _t0 = t0
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

    if do4029_resid:
        # for now, just lp12-502
        model_id = f"manual_20230617_mask_v0_nterms2"
        manual_csvpath = f'/Users/luke/Dropbox/proj/cpv/results/4029_mask/lc_lsresid_{model_id}.csv'
        df = pd.read_csv(manual_csvpath)
        times = np.array(df.time)
        # residual flux from subtracting the model given in model_id
        fluxs = np.array(df.r_flux)

    from astrobase.lcmath import find_lc_timegroups
    ngroups, groups = find_lc_timegroups(times, mingap=mingap)

    # Make plots
    plt.close('all')
    set_style("science")
    fig, ax = plt.subplots(figsize=(3,figsize_y))

    plot_period = np.nanmean(_periods)
    plot_t0 = t0s[0]
    plot_period_std = np.nanstd(_periods)

    ix = 0
    for n, g in enumerate(groups):

        cadence_sec = np.nanmedian(np.diff(times[g]))*24*60*60
        cadence_day = cadence_sec/(60*60*24)
        N_cycles_in_group = len(times[g]) * cadence_day / plot_period

        gtime = times[g]
        gflux = fluxs[g]

        # t = t0 + P*e
        e_start = int(np.floor((gtime[0] - plot_t0)/plot_period))
        e_end = int(np.floor((gtime[-1] - plot_t0)/plot_period))
        txt = f"{e_start} - {e_end}"

        if N_cycles_in_group <= 3:
            continue
        if len(gtime) < 100:
            continue
        #print(txt, N_cycles_in_group, len(gtime))

        if isinstance(cyclewindow, (tuple, list)):
            if e_start < cyclewindow[0]:
                continue
            if e_end > cyclewindow[1]:
                continue

        plot_phased_light_curve(
            gtime, gflux, plot_t0, plot_period, None,
            fig=fig, ax=ax,
            titlestr=None, binsize_phase=binsize_phase,
            xlim=xlim, yoffset=ix, showtext=txt, savethefigure=False,
            dy=yoffset, rasterized=rasterized
        )

        ix -= yoffset

    if do4029_resid:
        special_phases = [0.185, -0.29, -0.195, 0.265]
        for p in special_phases:
            ax.vlines(p, -100, 20, colors='darkgray', alpha=0.5,
                      linestyles='--', zorder=-10, linewidths=0.5)

    if showtitle:
        ax.set_title(
            f"{ticid.replace('_', ' ')} "
            f"P={plot_period*24:.4f}$\pm${plot_period_std*24:.3f}hr"
        )

    fig.text(-0.01,0.5, r"$\Delta$ Flux [%]", va='center',
             rotation=90)

    if isinstance(ylim, list):
        ax.set_ylim(ylim)

    ax.set_xlabel('Phase, φ')

    format_ax(ax)
    fig.tight_layout()

    s = ''
    if do4029_resid:
        s += 'resid_'
    if isinstance(cyclewindow, (tuple, list)):
        s += f'cycle{cyclewindow[0]}_to_{cyclewindow[1]}'
    _r = ''
    if rasterized:
        _r = '_rasterized'

    outpath = join(outdir, f"{s}{ticid}_P{plot_period*24:.4f}_{lc_cadences}_phase_timegroups{_r}.png")
    #fig.savefig(outpath, dpi=300)
    savefig(fig, outpath, dpi=600)


def plot_phase_timegroups_mosaic(
    outdir,
    ticid=None,
    lc_cadences='2min',
    manual_period=None,
    t0='binmin',
    ylim=None,
    binsize_phase=0.005,
    xlim=[-0.6,0.6],
    yoffset=5,
    showtitle=1,
    figsize_y=7,
    model_id=None,
    rasterized=False
    ):
    """
    As in plot_phase
    """

    assert ticid == "TIC_402980664", 'currently only tic4029 b/c of mosaic logic'

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
        if not isinstance(t0, float) and isinstance(manual_period, float):
            d = cpv_periodsearch(x_obs, y_flat, starid, outdir, t0=t0)
        else:
            d = {'times': x_obs, 'fluxs': y_flat,
                 't0': t0, 'period': manual_period}

        _t0 = d['t0']
        if isinstance(t0, float):
            _t0 = t0
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

    if isinstance(model_id, str):
        # for now, just lp12-502
        assert ticid == "TIC_402980664"
        manual_csvpath = f'/Users/luke/Dropbox/proj/cpv/results/4029_mask/lc_lsresid_{model_id}.csv'
        df = pd.read_csv(manual_csvpath)
        times = np.array(df.time)
        # residual flux from subtracting the model given in model_id
        fluxs = np.array(df.r_flux)

    from astrobase.lcmath import find_lc_timegroups
    ngroups, groups = find_lc_timegroups(times, mingap=3/24)

    # Make plots
    plt.close('all')
    set_style("science")
    factor=0.8
    factor=1.1
    fig, axs = plt.subplots(
        nrows=3, ncols=6,
        figsize=(factor*5.2, factor*(3.5/7)*7),
        constrained_layout=True
    )
    axs = axs.flatten()

    plot_period = np.nanmean(_periods)
    plot_t0 = t0s[0]
    plot_period_std = np.nanstd(_periods)

    ix = 0
    for n, g in enumerate(groups):

        ax = axs[ix]

        cadence_sec = np.nanmedian(np.diff(times[g]))*24*60*60
        cadence_day = cadence_sec/(60*60*24)
        N_cycles_in_group = len(times[g]) * cadence_day / plot_period

        gtime = times[g]
        gflux = fluxs[g]

        # t = t0 + P*e
        e_start = int(np.floor((gtime[0] - plot_t0)/plot_period))
        e_end = int(np.floor((gtime[-1] - plot_t0)/plot_period))
        txt = f"{e_start}"+"$\,$-$\,$"+f"{e_end}"

        if N_cycles_in_group <= 3:
            continue
        if len(gtime) < 100:
            continue
        #print(txt, N_cycles_in_group, len(gtime))

        plot_phased_light_curve(
            gtime, gflux, plot_t0, plot_period, None,
            fig=fig, ax=ax,
            binsize_phase=binsize_phase,
            xlim=xlim,
            #showtext=txt,
            titlestr=txt,
            titlepad=0.1,
            showtext=False,
            savethefigure=False,
            titlefontsize='x-small',
            rasterized=rasterized
        )
        if isinstance(ylim, (list, tuple)):
            ax.set_ylim(ylim)

        #ax.hlines(2.5, -0.6666/2, 0.6666/2, colors='darkgray', alpha=1,
        #          linestyles='-', zorder=-2, linewidths=1)

        ax.set_yticklabels([])
        ax.set_xticklabels([])

        gotyticks = False
        if isinstance(ylim, (list, tuple)):
            if max(ylim) >= 2 and min(ylim) <= -3:
                yticks = [-3,0,2]
                ax.set_yticks(yticks)
                gotyticks = True
        ax.minorticks_off()

        ax.set_xticks([-0.5, 0, 0.5])

        if ix in [0,6,12] and gotyticks:
            ax.set_yticklabels(yticks, fontsize='x-small')

        ix += 1

    for ix in range(12,18):
        axs[ix].set_xticklabels(['-0.5', '0', '0.5'], fontsize='x-small')

    axs[17].set_xticklabels([])

    if showtitle:
        fig.text(
            0.5, 1,
            f"LP 12-502 "
            f"(P={plot_period*24:.1f}h)",
            ha='center',
            va='center',
        )

    fig.text(-0.01,0.5, r"$\Delta$ Flux [%]", va='center', rotation=90, fontsize='large')
    fig.text(0.5,-0.01, r"Phase, φ", ha='center', fontsize='large')

    format_ax(ax)
    fig.tight_layout()
    axs[17].set_xticklabels(['-0.5', '0', '0.5'], fontsize='x-small')

    s = ''
    if isinstance(model_id, str):
        s += f'_{model_id}'
    if isinstance(ylim, (list, tuple)):
        s += f'_ymin{ylim[0]}_ymax{ylim[1]}'
    if rasterized:
        s += "_rasterized"

    outpath = join(
        outdir,
        f"{ticid}_P{plot_period*24:.4f}_{lc_cadences}_phase_timegroups_mosaic{s}.png"
    )

    fig.tight_layout(h_pad=0.2, w_pad=0.)

    fig.savefig(outpath, bbox_inches='tight', dpi=450)
    print(f"saved {outpath}")
    fig.savefig(outpath.replace('.png','.pdf'), bbox_inches='tight', dpi=600)
    print(f"saved {outpath.replace('.png','.pdf')}")


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

    fig.text(-0.01,0.5, "$\Delta$ Flux [%]", va='center',
             rotation=90)

    ax.set_xlabel('Time [TJD]')

    format_ax(ax)
    fig.tight_layout()

    savefig(fig, outpath, dpi=350)




def plot_phased_light_curve(
    time, flux, t0, period, outpath,
    ylim=None, xlim=[-0.6,0.6], binsize_phase=0.005, BINMS=2, titlestr=None,
    showtext=True, showtitle=False, figsize=None,
    c0='darkgray', alpha0=0.3,
    c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
    fig=None, ax=None, savethefigure=True,
    findpeaks_result=None, ylabel=None, showxticklabels=True,
    yoffset=0, dy=5, normfunc=True, xtxtoffset=0, titlepad=None,
    titlefontsize='small', showtimeticks=None, rasterized=False
    ):
    """
    Non-obvious args:
        binsize_phase (float): binsize in units of phase.
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
        showtimeticks: list of times to show with colored ticks
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

    x = time
    if normfunc:
        y = flux-np.nanmean(flux)
    else:
        y = flux

    # time units
    # x_fold = (x - t0 + 0.5 * period) % period - 0.5 * period
    # phase units
    # phase = (times - epoch)/period - floor((times - epoch)/period)
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
    if normfunc:
        norm = lambda x: 1e2*x + yoffset
    else:
        norm = lambda x: x

    if not plotnotscatter:
        ax.scatter(x_fold, norm(y), color=c0, marker='.',
                   s=1, rasterized=True, alpha=alpha0, linewidths=0)
    else:
        ax.plot(x_fold, norm(y), color=c0,
                lw=0.5, rasterized=True, alpha=alpha0)

    orb_bd = phase_bin_magseries(x_fold, y, binsize=binsize_phase, minbinelems=1)
    if c1 == 'k':
        ax.scatter(
            orb_bd['binnedphases'], norm(orb_bd['binnedmags']), color=c1,
            s=BINMS, linewidths=0,
            alpha=alpha1, zorder=1002, rasterized=rasterized#, linewidths=0.2, edgecolors='white'
        )
    else:
        ax.scatter(
            orb_bd['binnedphases'], norm(orb_bd['binnedmags']), color=c1,
            s=BINMS, linewidths=0.07, edgecolors='k',
            alpha=alpha1, zorder=1002, rasterized=rasterized#, linewidths=0.2, edgecolors='white'
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

    if isinstance(showtimeticks, (list, np.ndarray)):

        times = np.array(showtimeticks)
        phases = (times - t0)/period - np.floor((times - t0)/period)
        # wrap to [-0.5, 0.5]
        phases[phases > 0.5] -= 1

        if not xlim == [-0.6,0.6]:
            raise AssertionError("phasing off otherwise")

        tform = blended_transform_factory(ax.transData, ax.transAxes)

        ax.scatter(phases, np.ones(len(phases))*0.98+yoffset, transform=tform,
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
            if len(orb_bd['binnedmags'][sel]) > 0:
                ax.text(0.97+xtxtoffset,
                        np.nanmin(norm(orb_bd['binnedmags'][sel])-dy/5), txt,
                        transform=tform, ha='right',va='top', color='k',
                        fontsize=fontsize, bbox=props, zorder=999999)
            else:
                ax.text(0.97+xtxtoffset,
                        np.nanmin(norm(orb_bd['binnedmags'])-dy/5), txt,
                        transform=ax.transAxes, ha='right',va='top', color='k',
                        fontsize=fontsize, bbox=props, zorder=999999)

        else:
            if isinstance(t0, float):
                #txt = f'$t_0$ [BTJD]: {t0:.6f}\n$P$: {period:.6f} d'
                txt = f'{period*24:.1f} hr' # simpler label
                fontsize='small'
            elif isinstance(t0, int):
                txt = f'$t_0$ [BTJD]: {t0:.1f}\n$P$: {period:.6f} d'
                fontsize='xx-small'
            ax.text(0.97,0.03,txt,
                    transform=ax.transAxes,
                    ha='right',va='bottom', color='k', fontsize=fontsize)

    if showtitle:
        txt = f'$t_0$ [BTJD]: {t0:.6f}. $P$: {period:.6f} d'
        ax.set_title(txt, fontsize='small', pad=titlepad)

    if isinstance(titlestr,str):
        ax.set_title(titlestr.replace("_"," "), fontsize=titlefontsize, pad=titlepad)

    if savethefigure:
        ax.set_ylabel("$\Delta$ Flux [%]")
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
    else:
        if showxticklabels == True:
            pass
        elif isinstance(showxticklabels, (list, tuple)):
            ax.set_xticks(showxticklabels)
            ax.set_xticklabels(showxticklabels)

    if isinstance(ylabel, str):
        ax.set_ylabel(ylabel)

    if fig is not None:
        fig.tight_layout()

    if savethefigure:
        ax.spines[['right', 'top']].set_visible(False)
        savefig(fig, outpath, dpi=350)
        plt.close('all')
    else:
        if fig is not None:
            return fig, ax
        else:
            pass


def get_ylimguess(y):
    ylow = np.nanpercentile(y, 2.5)
    yhigh = np.nanpercentile(y, 97.5)
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
    binsize_phase=0.005,
    lcpipeline='spoc2min'
    ):

    # get data
    d = periodsearch_result
    r = findpeaks_result
    hdul = fits.open(lcpath)
    hdr = hdul[0].header
    data = hdul[1].data
    hdul.close()

    # quality flags
    QUALITYKEYDICT = {
        'spoc2min': 'QUALITY',
        'qlp': 'QUALITY',
        'cdips': 'IRQ3'
    }
    qual = data[QUALITYKEYDICT[lcpipeline]]
    if lcpipeline in ['spoc2min', 'qlp']:
        sel = (qual == 0)
    elif lcpipeline == 'cdips':
        sel = (qual == 'G')

    # centroid data
    CENTRKEYDICT = {
        'spoc2min': [
            'MOM_CENTR2', # column
            'MOM_CENTR1' # row
        ],
        'qlp': [
            'SAP_X', # column
            'SAP_Y' # row
        ],
        'cdips': [
            'XIC', # column
            'YIC' # row
        ]
    }
    xc = data[CENTRKEYDICT[lcpipeline][0]][sel] # column
    yc = data[CENTRKEYDICT[lcpipeline][1]][sel] # row

    # background data
    BKGDKEYDICT = {
        'spoc2min': 'SAP_BKG',
        'qlp': 'SAP_BKG',
        'cdips': 'BGV'
    }
    bgv = data[BKGDKEYDICT[lcpipeline]][sel]

    assert len(xc) == len(d['times'])

    # make plot
    plt.close('all')
    set_style("clean")

    #fig = plt.figure(figsize=(8,4.5))
    fig = plt.figure(figsize=(8,6))
    axd = fig.subplot_mosaic(
        """
        AAAABBCC
        AAAABBCC
        AAAADDDD
        EEFFKKLL
        GGHHKKLL
        IIJJMMLL
        """
    )

    axd['E'].get_shared_x_axes().join(axd['E'], axd['G'])
    axd['E'].get_shared_x_axes().join(axd['E'], axd['I'])
    axd['G'].get_shared_x_axes().join(axd['G'], axd['I'])

    axd['F'].get_shared_x_axes().join(axd['F'], axd['H'])
    axd['F'].get_shared_x_axes().join(axd['F'], axd['J'])
    axd['H'].get_shared_x_axes().join(axd['H'], axd['J'])


    # pdcsap flux vs time (qual==0, after 5-day median smooth)
    ax = axd['D']
    bd = time_bin_magseries(d['times'], d['fluxs'], binsize=1200, minbinelems=1)
    #ax.scatter(d['times'], d['fluxs'], c='lightgray', s=1, zorder=1)
    ax.scatter(bd['binnedtimes'], bd['binnedmags'], c='k', s=0.2, zorder=2,
               rasterized=True)
    #txt = 'PDCSAP, 5d median filter'
    #bbox = dict(facecolor='white', alpha=0.9, pad=0, edgecolor='white')
    #ax.text(0.03, 0.97, txt, ha='left', va='top', bbox=bbox, zorder=3,
    #        transform=ax.transAxes)
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
    P_harmonics = []
    for ix in range(1,11):
        P_harmonics.append(ix*d['period'])
        P_harmonics.append(d['period']/ix)

    sel = (
        (P_harmonics > np.nanmin(d['lsp']['periods']))
        &
        (P_harmonics < np.nanmax(d['lsp']['periods']))
    )
    ax.vlines(nparr(P_harmonics)[sel], ymin, ymax, colors='darkgray',
              alpha=0.5, linestyles=':', zorder=-2, linewidths=0.5)
    ax.set_ylim([ymin, ymax])
    ax.update({'xlabel': 'Period [d]', 'ylabel': 'PDM Θ', 'xscale': 'log'})

    # phased LC
    ax = axd['A']
    ylim = get_ylimguess(1e2*(bd['binnedmags']-np.nanmean(bd['binnedmags'])))
    plot_phased_light_curve(
        d['times'], d['fluxs'], d['t0'], d['period'], None, ylim=ylim,
        xlim=[-0.6,0.6], binsize_phase=0.005, BINMS=2, titlestr=None,
        showtext=True, showtitle=False, figsize=None, c0='darkgray',
        alpha0=0.3, c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
        fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
        showxticklabels=True
    )
    ax.set_ylabel("$\Delta$ Flux [%]")
    ax.set_xlabel("Phase, φ")

    # phased LC at 2x period
    ax = axd['G']
    plot_phased_light_curve(
        d['times'], d['fluxs'], d['t0'], 2*d['period'], None, ylim=ylim,
        xlim=[-0.6,0.6], binsize_phase=0.005, BINMS=2, titlestr=None,
        showtext=False, showtitle=False, figsize=None, c0='darkgray',
        alpha0=0.3, c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
        fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
        showxticklabels=True, ylabel=r'2$\times P$'
    )
    ax.set_xticklabels(['-0.5', '', '0', '', '0.5'])
    ax.set_xlabel("φ")

    # phased LC at 0.5x period
    ax = axd['I']
    plot_phased_light_curve(
        d['times'], d['fluxs'], d['t0'], 0.5*d['period'], None, ylim=ylim,
        xlim=[-0.6,0.6], binsize_phase=0.005, BINMS=2, titlestr=None,
        showtext=False, showtitle=False, figsize=None, c0='darkgray',
        alpha0=0.3, c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
        fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
        showxticklabels=True, ylabel=r'0.5$\times P$'
    )
    ax.set_xticklabels(['-0.5', '', '0', '', '0.5'])
    ax.set_xlabel("φ")

    # # phased norm LC, w/ smooth model
    ax = axd['E']
    norm = lambda x: 1e2 * (x - np.nanmean(x))
    ax.scatter(
        r['binned_phase'], norm(r['binned_orig_flux']), c='k', zorder=1, s=1.5
    )
    ax.plot(
        r['binned_phase'], norm(r['binned_trend_flux']), c='darkgray', zorder=2
    )
    tform = blended_transform_factory(ax.transData, ax.transAxes)
    peakloc = findpeaks_result['peaks_phaseunits']
    peakloc[peakloc>0.5] -= 1 # there must be a smarter way
    ax.scatter(peakloc, np.ones(len(peakloc))*0.98, transform=tform,
               marker='v', s=10, linewidths=0, edgecolors='none',
               color='k', alpha=0.5, zorder=1000)
    ax.set_xlim([-0.6,0.6])
    ax.set_xticks([-0.5, -0.25, 0, 0.25, 0.5])
    ax.set_xticklabels(['-0.5', '', '0', '', '0.5'])
    ax.set_xlabel("φ")
    ax.set_ylabel('dips')


    # # resids of above, w/ ticks for auto dip find
    # ax = axd['G']

    # # phased flux in BGD aperture
    ax = axd['F']
    nbgv = bgv/np.nanmedian(bgv)
    sel_bgv = np.isfinite(nbgv) & (nbgv < np.nanpercentile(nbgv, 98))
    ylim = get_ylimguess(nbgv[sel_bgv])
    ylim = [0.9,1.1]
    try:
        plot_phased_light_curve(
            d['times'][sel_bgv], nbgv[sel_bgv], d['t0'], d['period'], None,
            ylim=ylim,
            xlim=[-0.6,0.6], binsize_phase=0.01, BINMS=2, titlestr=None,
            showtext=False, showtitle=False, figsize=None, c0='darkgray',
            alpha0=0.3, c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
            fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
            showxticklabels=True, ylabel='Bkg', normfunc=False
        )
    except TypeError:
        pass
    ax.set_xticklabels(['-0.5', '', '0', '', '0.5'])
    ax.set_xlabel("φ")

    #
    # # phased ctd x/y
    #
    ax = axd['H']
    if np.any(np.isfinite(xc)):
        plot_phased_light_curve(
            d['times'], xc, d['t0'], d['period'], None, ylim=None,
            xlim=[-0.6,0.6], binsize_phase=0.01, BINMS=2, titlestr=None,
            showtext=False, showtitle=False, figsize=None, c0='C0',
            alpha0=0.2, c1='C0', alpha1=1, phasewrap=True, plotnotscatter=False,
            fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
            showxticklabels=True, ylabel='c$_\mathrm{col}$'
        )
    ax.set_xticklabels(['-0.5', '', '0', '', '0.5'])
    ax.set_xlabel("φ")
    ax = axd['J']
    if np.any(np.isfinite(yc)):
        plot_phased_light_curve(
            d['times'], yc, d['t0'], d['period'], None, ylim=None,
            xlim=[-0.6,0.6], binsize_phase=0.01, BINMS=2, titlestr=None,
            showtext=False, showtitle=False, figsize=None, c0='C1',
            alpha0=0.2, c1='C1', alpha1=1, phasewrap=True, plotnotscatter=False,
            fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
            showxticklabels=True, ylabel='c$_\mathrm{row}$'
        )
    ax.set_xticklabels(['-0.5', '', '0', '', '0.5'])
    ax.set_xlabel("φ")

    # join subplots

    #
    # # star info (Gaia, TIC8, dip search)
    #
    ax = axd['L']
    ax.set_axis_off()

    # tic8 info
    TEFFKEYDICT = {
        'spoc2min': 'TEFF',
        'qlp': 'TEFF',
        'cdips': 'TICTEFF'
    }
    PMRAKEYDICT = {
        'spoc2min': 'PMRA',
        'qlp': 'PMRA',
        'cdips': 'PM_RA[mas/yr]'
    }
    PMDECKEYDICT = {
        'spoc2min': 'PMDEC',
        'qlp': 'PMDEC',
        'cdips': 'PM_Dec[mas/year]'
    }

    ticid = str(hdr['TICID'])
    sector = str(hdr['SECTOR'])
    cam = str(hdr['CAMERA'])
    ccd = str(hdr['CCD'])
    Tmag = f"{hdr['TESSMAG']:.1f}"
    if hdr[TEFFKEYDICT[lcpipeline]] is not None and hdr[TEFFKEYDICT[lcpipeline]] != 'nan':
        teff_tic8 = f"{int(hdr[TEFFKEYDICT[lcpipeline]]):d} K"
    else:
        teff_tic8 = f"NaN"
    if hdr['RA_OBJ'] != 'nan':
        ra = f"{hdr['RA_OBJ']:.2f}"
    else:
        ra = 'NaN'
    if hdr['DEC_OBJ'] != 'nan':
        dec = f"{hdr['DEC_OBJ']:.2f}"
    else:
        dec = 'NaN'
    if hdr[PMRAKEYDICT[lcpipeline]] is not None and hdr[PMRAKEYDICT[lcpipeline]] != 'nan':
        pmra = f"{hdr[PMRAKEYDICT[lcpipeline]]:.1f}"
    else:
        pmra = 'NaN'
    if hdr[PMDECKEYDICT[lcpipeline]] is not None and hdr[PMDECKEYDICT[lcpipeline]] != 'nan':
        pmdec = f"{hdr[PMDECKEYDICT[lcpipeline]]:.1f}"
    else:
        pmdec = 'NaN'
    ra_obj, dec_obj = hdr['RA_OBJ'], hdr['DEC_OBJ']
    c_obj = SkyCoord(ra_obj, dec_obj, unit=(u.deg), frame='icrs')

    # gaia info
    dr2_source_id = tic_to_gaiadr2(ticid)

    runid = f"dr2_{dr2_source_id}"

    dr2_source_ids = np.array([np.int64(dr2_source_id)])
    try:
        gdf = given_source_ids_get_gaia_data(
            dr2_source_ids, runid, n_max=5, overwrite=False,
            enforce_all_sourceids_viable=True, savstr='', which_columns='*',
            table_name='gaia_source', gaia_datarelease='gaiadr2', getdr2ruwe=False
        )
    except Exception as e:
        LOGERROR(f'{ticid} = {runid} failed due to:\n{e}...')
        return 1
    try:
        gdf_ruwe = given_source_ids_get_gaia_data(
            dr2_source_ids, runid+"_ruwe", n_max=5, overwrite=False,
            enforce_all_sourceids_viable=True, savstr='', which_columns='*',
            table_name='gaia_source', gaia_datarelease='gaiadr2', getdr2ruwe=True
        )
        ruwe =  f"{gdf_ruwe['ruwe'].iloc[0]:.2f}"
    except AttributeError:
        ruwe = 'NaN'

    Gmag = f"{gdf['phot_g_mean_mag'].iloc[0]:.1f}"
    Rpmag = f"{gdf['phot_rp_mean_mag'].iloc[0]:.1f}"
    Bpmag = f"{gdf['phot_bp_mean_mag'].iloc[0]:.1f}"
    bp_rp = f"{gdf['bp_rp'].iloc[0]:.2f}"
    plx = f"{gdf['parallax'].iloc[0]:.2f}"
    plx_err = f"{gdf['parallax_error'].iloc[0]:.2f}"
    dist_pc = 1/(float(plx)*1e-3)
    dist = f"{dist_pc:.1f}"

    # nbhr info
    ticids, tmags = get_2px_neighbors(c_obj, hdr['TESSMAG'])
    brightest_inds_first = np.argsort(tmags)
    ticids = ticids[brightest_inds_first]
    tmags = tmags[brightest_inds_first]
    N_nbhrs = len(ticids)-1
    nbhrstr = ''
    MAX_N = 7
    if N_nbhrs >= 1:
        for _ticid, tmag in zip(ticids[1:MAX_N], tmags[1:MAX_N]):
            nbhrstr += f"TIC {_ticid}: ΔT={tmag-hdr['TESSMAG']:.1f}\n"

    txt = (
        f"TIC {ticid}\n"
        f"GDR2 {dr2_source_id}\n"
        f"SEC{sector}, CAM{cam}, CCD{ccd}\n"
        "—\n"
        f"α={ra}, δ={dec} (deg)\n"
        r"$\mu_\alpha$="+pmra+r", $\mu_\delta$="+pmdec+f" (mas/yr)\n"
        f"T={Tmag}\n"
        f"TIC8 Teff={teff_tic8}\n"
        f"G={Gmag}, RP={Rpmag}, BP={Bpmag}\n"
        f"BP-RP={bp_rp}\n"
        f"RUWE={ruwe}\n"
        f"plx={plx}"+"$\pm$"+f"{plx_err} mas\n"
        f"d={dist} pc\n"
        "—\n"
        'N$_{\mathrm{dip}}$='+f'{r["N_peaks"]}\n'
        f'P2P={1e2*r["p2p_est"]:.1f}%\n'
        r'δ$_{\mathrm{thresh}}$'+f'={1e2*r["height"]:.1f}%\n'
        f'φ={r["peaks_phaseunits"]}\n'
        "—\n"
        'N$_{\mathrm{nbhr}}$: '+f'{N_nbhrs}\n'
        f'{nbhrstr}\n'
    )

    txt_x = -0.1
    txt_y = 0.5
    ax.text(txt_x, txt_y, txt, ha='left', va='center', fontsize='small', zorder=2,
            transform=ax.transAxes)

    #
    # gaia CAMD
    #
    ss = axd['K'].get_subplotspec()
    axd['K'].remove()
    import mpl_scatter_density # adds projection='scatter_density'
    axd['K'] = fig.add_subplot(ss, projection='scatter_density')
    ax = axd['K']

    get_xval_no_corr = (
        lambda _df: np.array(_df['phot_bp_mean_mag'] - _df['phot_rp_mean_mag'])
    )
    get_yval_no_corr = (
        lambda _df: np.array(
            _df['phot_g_mean_mag'] + 5*np.log10(_df['parallax']/1e3) + 5
        )
    )

    underplot_gcns(ax, get_xval_no_corr, get_yval_no_corr)

    bp_rp = get_xval_no_corr(gdf)
    M_G = get_yval_no_corr(gdf)
    ax.scatter(bp_rp, M_G, zorder=9000,
               c='lime', s=60, marker='X', edgecolors='k', linewidths=0.5)

    xmin, xmax = 0.9, 3.6
    if bp_rp > xmax:
        xmax = bp_rp + 0.2
    ymax, ymin = 11.5, 5
    if M_G > ymax:
        ymax = M_G + 0.5
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymax, ymin])
    ax.set_xlabel('BP-RP')
    ax.set_ylabel(r'M$_{\rm G}$')

    #
    # banyan overlay
    #
    ax = axd["M"]
    ax.set_axis_off()

    import sys
    from complexrotators.paths import BANYANDIR
    sys.path.append(BANYANDIR)
    from core import membership_probability

    ra, dec = float(gdf.ra), float(gdf.dec)
    pmra, pmdec = float(gdf.pmra), float(gdf.pmdec)
    epmra, epmdec = float(gdf.pmra_error), float(gdf.pmdec_error)
    plx, eplx = float(gdf.parallax), float(gdf.parallax_error)
    output = membership_probability(ra=ra, dec=dec, pmra=pmra, pmdec=pmdec,
                                    epmra=epmra, epmdec=epmdec, plx=plx, eplx=eplx,
                                    use_plx=True, use_rv=False)

    probs = np.array(output['ALL'].iloc[0].round(4))
    assocs = np.array(output['ALL'].iloc[0].index)
    banyan_df = pd.DataFrame({"prob": probs}, index=assocs)

    sel = banyan_df.prob > 1e-4
    sdf = banyan_df[sel].sort_values(by='prob', ascending=False)

    csvdir = os.path.dirname(outpath)
    csvname = f"TIC{ticid}_GDR2{dr2_source_id}_banyan_result.csv"
    csvpath = join(csvdir, csvname)
    sdf.to_csv(csvpath)

    txt = sdf.__repr__()
    txt_x = 0.5
    txt_y = 0.5
    ax.text(txt_x, txt_y, txt, ha='center', va='center', fontsize='small', zorder=2,
            transform=ax.transAxes)

    #
    # DSS query
    #
    ra = hdr['RA_OBJ']
    dec = hdr['DEC_OBJ']

    dss_overlay(fig, axd, ra, dec)

    # set naming options
    s = ''

    # height / weight
    fig.tight_layout(h_pad=0)

    savefig(fig, outpath)


def dss_overlay(fig, axd, ra, dec):

    from astrobase.plotbase import skyview_stamp
    from astropy.wcs import WCS

    sizepix = 220
    try:
        dss, dss_hdr = skyview_stamp(ra, dec, survey='DSS2 Red',
                                     scaling='Linear', convolvewith=None,
                                     sizepix=sizepix, flip=False,
                                     cachedir='~/.astrobase/stamp-cache',
                                     verbose=True, savewcsheader=True)
    except (OSError, IndexError, TypeError) as e:
        LOGERROR('downloaded FITS appears to be corrupt, retrying...')
        try:
            dss, dss_hdr = skyview_stamp(ra, dec, survey='DSS2 Red',
                                         scaling='Linear', convolvewith=None,
                                         sizepix=sizepix, flip=False,
                                         cachedir='~/.astrobase/stamp-cache',
                                         verbose=True, savewcsheader=True,
                                         forcefetch=True)

        except Exception as e:
            LOGERROR('failed to get DSS stamp ra {} dec {}, error was {}'.
                     format(ra, dec, repr(e)))
            return None, None

    ss = axd['C'].get_subplotspec()
    axd['C'].remove()
    axd['C'] = fig.add_subplot(ss, projection=WCS(dss_hdr))
    ax = axd['C']

    cset = ax.imshow(dss, origin='lower', cmap=plt.cm.gray_r, zorder=-2)

    ax.set_xlabel(' ')
    ax.set_ylabel(' ')

    ax.grid(ls='--', alpha=0.5)

    # DSS is ~1 arcsecond per pixel. overplot 1px and 2px apertures
    for ix, radius_px in enumerate([21, 21*2]):
        circle = plt.Circle(
            (sizepix/2, sizepix/2), radius_px, color='C{}'.format(ix),
            fill=False, zorder=5+ix, lw=0.5, alpha=0.5
        )
        ax.add_artist(circle)

    props = dict(boxstyle='square', facecolor='lightgray', alpha=0.3, pad=0.15,
                 linewidth=0)
    ax.text(0.97, 0.03, 'r={1,2}px', transform=ax.transAxes,
            ha='right',va='bottom', color='k', bbox=props, fontsize='x-small')



def get_2px_neighbors(c_obj, tessmag):

    radius = 42*u.arcsecond

    #nbhr_stars = Catalogs.query_region(
    #    "{} {}".format(float(c_obj.ra.value), float(c_obj.dec.value)),
    #    catalog="TIC",
    #    radius=radius
    #)
    nbhr_stars = Catalogs.query_region(
        c_obj,
        catalog="TIC",
        radius=radius
    )

    ticids = nbhr_stars[nbhr_stars['Tmag'] < tessmag+2.5]['ID']
    tmags = nbhr_stars[nbhr_stars['Tmag'] < tessmag+2.5]['Tmag']

    return ticids, tmags



def underplot_gcns(ax, get_xval_no_corr, get_yval_no_corr):

    from rudolf.helpers import get_gaia_catalog_of_nearby_stars
    df_bkgd = get_gaia_catalog_of_nearby_stars()

    from matplotlib.colors import LinearSegmentedColormap
    # "Viridis-like" colormap with white background
    white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
        (0, '#ffffff'),
        (1e-20, '#440053'),
        (0.2, '#404388'),
        (0.4, '#2a788e'),
        (0.6, '#21a784'),
        (0.8, '#78d151'),
        (1, '#fde624'),
    ], N=256)

    _x = get_xval_no_corr(df_bkgd)
    _y = get_yval_no_corr(df_bkgd)
    s = np.isfinite(_x) & np.isfinite(_y)

    # add log stretch...
    from astropy.visualization import LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    vmin, vmax = 10, 1000

    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())

    cmap = "Greys"
    #cmap = white_viridis
    density = ax.scatter_density(_x[s], _y[s], cmap=cmap, norm=norm, zorder=1)


def get_flx_wav_given_2d_and_target(flx_2d, wav_2d, target_wav):
    _preorder = np.argmin(np.abs(wav_2d - target_wav), axis=1)
    viable_orders = np.argwhere(
        (_preorder != wav_2d.shape[1]-1) & (_preorder != 0)
    )
    order = int(
        viable_orders[np.argmin(
            np.abs(_preorder[viable_orders] - wav_2d.shape[1]/2)
        )]
    )

    flx, wav = flx_2d[order, :], wav_2d[order, :]
    return flx, wav


def plot_spectrum_windows(outdir, starid, inst='HIRES', ylims=None):

    # note: could add in K 7699..
    #lines = ['Ca K', 'Ca H & Hε', 'He', 'Hα', 'Li', 'K']
    lines = ['Ca K', 'Ca H & Hε', 'Hα', 'Li']

    deltawav = 7.5
    xlims = [
        [3933.66-deltawav, 3933.66+deltawav], # ca k
        [3968.47-deltawav, 3968.47+deltawav], # ca h
        #[5875.618-deltawav, 5875.618+deltawav], # He
        #[5895.92-1.5*deltawav, 5895.92+deltawav], # Na D1
        [6562.8-deltawav, 6562.8+deltawav], # halpha
        [6707.8-deltawav, 6707.8+deltawav], # li6708
        #[7699.-deltawav, 7699+deltawav], # K
    ]
    if ylims is None:
        ylims = [
            None,#[-1, 35],
            None,#[-1, 35],
            #[0.8, 1.95],
            None,
            None,
            #None,
        ]
    xticks = [
        [3930, 3940],
        [3965, 3975],
        [6560, 6570],
        #None,
        [6705, 6715],
        #[7695, 7705],
    ]
    globs = [
        '*bj*',#order07*', # ca k
        '*bj*',#order07*', # ca h
        #'*rj*',#order10*', # he
        #'*rj*',#order11*', # Na D1
        '*ij*',#order00*', # H alpha
        '*ij*',#order01*', # Li6708
        #'*ij*',#orderXX*', # K
    ]

    #
    # make plot
    #
    plt.close('all')
    set_style('science')

    f = 1.5
    fig, axs = plt.subplots(ncols=len(lines), figsize=(f*3.5,f*1.15))
    #fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(2.5,2.5))
    axs = axs.flatten()

    from cdips_followup.spectools import read_hires

    for ix, line, xlim, ylim, xtick, _glob in zip(
        range(len(lines)), lines, xlims, ylims, xticks, globs
    ):

        fitspaths = glob(join(DATADIR, 'spectra', inst, starid, _glob))
        assert len(fitspaths) == 1
        spectrum_path = fitspaths[0]

        hdul = fits.open(spectrum_path)
        dateobs = hdul[0].header['DATE-OBS']
        exptime = hdul[0].header['EXPTIME']
        utc = hdul[0].header['UTC'][:5]
        timestr = f"{dateobs} {utc}UT"
        hdul.close()

        flx_2d, wav_2d = read_hires(
            spectrum_path, is_registered=0, return_err=0
        )

        ax = axs[ix]

        norm = lambda x: x/np.nanmedian(x)
        fn = lambda x: gaussian_filter1d(x, sigma=2)

        flx, wav = get_flx_wav_given_2d_and_target(
            flx_2d, wav_2d, xlim[1]-deltawav
        )

        sel = ((xlim[0]-20) < wav) & (wav < xlim[1]+20)
        xv = wav[sel]
        yv = fn(norm(flx[sel]))
        ax.plot(
            xv, yv, c='k', zorder=3, lw=0.2
        )

        sel0 = ((xlim[0]-0) < wav) & (wav < xlim[1]+0)
        max_y = 1.5*np.nanmax(fn(norm(flx[sel0])))

        ax.set_title(line, fontsize='small', pad=0)
        ax.set_xlim(xlim)

        if ylim is not None:
            ax.set_ylim(ylim)
        else:
            ax.set_ylim([0, max_y])

        if xtick is not None:
            ax.set_xticks(xtick)
            ax.set_xticklabels([str(x) for x in xtick])

        if line == 'K':
            ax.yaxis.set_major_locator(MultipleLocator(0.5))
        if line == 'Li':
            ax.yaxis.set_major_locator(MultipleLocator(0.5))
        if line == 'Hα':
            ax.yaxis.set_major_locator(MultipleLocator(2))



    #axs[0].set_ylabel("Relative flux")
    stitle = f'{starid.replace("TIC","TIC ")} {timestr}'
    print(stitle + ' ' + str(exptime))
    fig.text(0.5, 1.01, stitle, va='top', ha='center', rotation=0)
    fig.text(-0.01,0.5, r'Relative flux', va='center', ha='center',
             rotation=90)
    fig.text(0.5,-0.01, r'Wavelength [$\AA$]', va='center', ha='center',
             rotation=0)

    #fig.tight_layout(w_pad=0.2, h_pad=0.5)
    fig.tight_layout(w_pad=0.2)

    # set naming options
    s = ''

    outpath = os.path.join(outdir, f'{starid}_{inst}_spectrum_windows{s}.png')
    savefig(fig, outpath, dpi=400)


def plot_winered_windows(outdir, inst='WINERED',
                         datestr='june10', frameid='WINA00038986',
                         starid='TIC_89026133'):
    """
    atomic Ti (m=167, 168, 180, 182), and Ca (m=171)

    m165 is REALLY GOOD!!!
        [10661, 10669]

    m167 needs decent sky subtrxn

    m168 is great, no contamination at all
        [10496, 10504]

    m171 Ca is nice, but beware extraction near chip edges
        [10344, 10350.5]

    m181 is a bit crappy -- skip

    m182 has two nice lines, tho 9708 needs sky sub
        [9676, 9682]
        [9706, 9712]

    """

    if datestr == 'june10':
        vnum = 'v1'
    elif datestr == 'june03':
        vnum = 'v0'

    lines = 'Ti Ti Ca Ti Ti'.split()

    wavs = [
        9678.346822347512, # Ti m182 first
        9708.478432545879, # Ti m182 second
        10346.80091026195, #  Ca m171
        10499.156483768491, # Ti m168
        10664.70983956067, # Ti m165
    ]

    deltawav = 5 # angstrom
    xlims = [
        [w-deltawav, w+deltawav] for w in wavs
    ]
    ylims = [
        [0.8,1.1],
        [0.8,1.1],
        [0.8,1.1],
        [0.8,1.1],
        [0.75,1.1]
    ]
    xticks = [
        [9675, 9680],
        [9705, 9710],
        [10345, 10350],
        [10496, 10501],
        [10661, 10666]
    ]
    globs = [
        '*m182*',
        '*m182*',
        '*m171*',
        '*m168*',
        '*m165*',
    ]

    badbands = [
        [9707.9,9708.3],
        [9711.2,9711.6],
        [9703.5,9705.],
        [10349.7,10350.6]
    ]

    #
    # make plot
    #
    plt.close('all')
    #set_style('clean')
    set_style('science')

    f = 2.5
    fig, axs = plt.subplots(ncols=len(lines), figsize=(f*3.3,f*1.15))
    axs = axs.flatten()

    from cdips_followup.spectools import read_winered

    def get_lam_to_dv(wav0):

        def lam_to_dv(lam):
            dv_by_c = (lam - wav0)/wav0
            dv = dv_by_c * 2.998e5 # units: km/s
            return dv

        return lam_to_dv

    def get_dv_to_lam(wav0):

        def dv_to_lam(dv):
            dv_by_c = dv / 2.998e5
            lam = (dv_by_c * wav0) + wav0
            return lam

        return dv_to_lam


    for ix, line, xlim, ylim, xtick, _glob in zip(
        range(len(lines)), lines, xlims, ylims, xticks, globs
    ):

        airvac = "VAC"
        fitspaths = glob(join(DATADIR, 'spectra', inst,
                              f"bouma_{datestr}_hiresy",
                              f"{starid}_{frameid}_output_{vnum}", f"{starid}_sum",
                              f"{airvac}_NORM", "fsr1.05", _glob))
        assert len(fitspaths) == 1
        spectrum_path = fitspaths[0]

        flx, wav = read_winered(spectrum_path)

        ax = axs[ix]

        norm = lambda x: x/np.nanmedian(x)
        fn = lambda x: gaussian_filter1d(x, sigma=4)

        sel = ((xlim[0]-20) < wav) & (wav < xlim[1]+20)
        ax.plot(
            wav[sel], fn(norm(flx[sel])), c='k', zorder=3, lw=0.2
        )

        ax.set_title(line)
        ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        if xtick is not None:
            ax.set_xticks(xtick)
            ax.set_xticklabels([str(x) for x in xtick])

        for badband in badbands:
            ymin, ymax = ax.get_ylim()
            _y = np.linspace(ymin, ymax, 100)
            ax.fill_betweenx(_y, np.ones_like(_y)*badband[0],
                             np.ones_like(_y)*badband[1], zorder=-100, color='k',
                             alpha=0.2, ec='none')

        wav = wavs[ix]
        ax.vlines(wav, ymin, ymax, colors='k',
                  linestyles=':', zorder=-5, lw=0.5, alpha=0.5)

        ax.set_ylim(ylim)

    for ix in range(len(lines)):
        wav0 = wavs[ix]
        secax_x = axs[ix].secondary_xaxis(
            -0.2, functions=(get_lam_to_dv(wav0), get_dv_to_lam(wav0))
        )



    fig.text(-0.01,0.5, r'Relative flux', va='center', ha='center',
             rotation=90)

    fig.text(0.5,-0.01, r'Wavelength [$\AA$] and Velocity [km/s]', va='center', ha='center',
             rotation=0)

    fig.tight_layout(w_pad=0.2)

    # set naming options
    s = ''

    outpath = os.path.join(outdir, f'{starid}_{inst}_{frameid}_spectrum_windows{s}.png')
    savefig(fig, outpath, dpi=400)



def plot_quasiperiodic_removal_diagnostic(d, pngpath):

    gfn = lambda x: gaussian_filter1d(x, sigma=1)

    # d = {
    #    'best_interp_key': best_interp_key,
    #    'pgdict': pgdict, # nominal periodogram dictionary
    #    'resid_pgdict': resid_pgdict, # periodogram on the RESIDUAL dictionary
    #    'period': period,
    #    'fn': fn,
    #    'time': time, # input time
    #    'flux': flux, # input flux
    #    'x_sw': x_sw, # sorted & wrapped phase
    #    'y_sw': y_sw, # sorted & wrapped flux
    #    'x_nsnw': x_nsnw, # not sorted & not wrapped phase
    #    'y_nsnw': y_nsnw, # not sorted & not wrapped flux
    #    'y_model_sw': y_model_sw, # sorted & wrapped model flux
    #    'y_model_nsnw': y_model_nsnw, # not sorted & not wrapped model flux
    #    'y_resid_sw': y_resid_sw, # y_sw - y_model_sw
    #    'y_resid_nsnw': y_resid_nsnw, # y_nsnw - y_model_nsnw
    #    'x_resid_sw_b': x_resid_sw_b,
    #    'y_resid_sw_b': y_resid_sw_b,
    #    'x_resid_nsnw_b': x_resid_nsnw_b,
    #    'y_resid_nsnw_b': y_resid_nsnw_b,
    #    'x_w_model': x_w_model,
    #    'y_w_model': y_w_model,
    #    'x_nw_model': x_nw_model,
    #    'y_nw_model': y_nw_model,
    #}

    # NOTE: populate the namespace using all keys in the passed dictionary.
    from types import SimpleNamespace
    ns = SimpleNamespace(**d)
    pgdict = ns.pgdict
    resid_pgdict = ns.resid_pgdict

    plt.close('all')
    set_style('clean')

    fig = plt.figure(figsize=(9,6), layout='constrained')
    axd = fig.subplot_mosaic(
        #"""
        #AABBEE
        #AABBEE
        #CCCCCC
        #DDDDDD
        #"""
        """
        AAACCCCCC
        AAACCCCCC
        BBBCCCCCC
        BBBDDDDDD
        EEEDDDDDD
        EEEDDDDDD
        FFFGGGHHH
        FFFGGGHHH
        IIIJJJKKK
        IIIJJJKKK
        """
    )

    n = lambda x: 1e2 * (x - np.nanmean(x))

    # flux vs phase
    ax = axd['A']
    ax.scatter(ns.x_sw, n(ns.y_sw), zorder=1, s=0.2, c='lightgray', linewidths=0)
    ax.scatter(ns.x_sw_b, n(ns.y_sw_b), zorder=3, s=1, c='k', linewidths=0)
    ax.plot(ns.x_w_model, n(ns.y_w_model), lw=0.5, zorder=2, c='C0')
    ylim = get_ylimguess(n(ns.y_sw_b))
    ax.set_ylim(ylim)
    ax.update({#'xlabel': 'φ',
               'ylabel': '$\Delta$ Flux [%]', 'xlim':[-0.6,0.6],
               'title': f"P = {pgdict['period']*24:.2f} hr"
              })
    ax.set_xticklabels([])

    # resid flux vs phase
    ax = axd['B']
    ax.scatter(ns.x_sw, n(ns.y_resid_sw), zorder=1, s=0.2, c='lightgray',
               linewidths=0)
    ax.scatter(ns.x_resid_sw_b, n(ns.y_resid_sw_b), zorder=3, s=1, c='k',
               linewidths=0)
    ax.plot(ns.x_w_model, n(ns.y_w_model-ns.y_w_model), lw=0.5, zorder=2, c='C0')
    factor = 15
    ylim = factor*np.array(get_ylimguess(n(ns.y_resid_sw_b)))
    ax.set_ylim(ylim)
    ax.update({#'xlabel': 'φ',
               'ylabel': 'Resid [%]', 'xlim':[-0.6,0.6],
              })
    ax.set_xticklabels([])

    # resid flux vs phase colored and binned by time
    ax = axd['E']
    cmap = mpl.colormaps['Spectral']
    norm = mpl.colors.Normalize(vmin=ns.cyclenum_sw.min(), vmax=ns.cyclenum_sw.max())
    binsize_minutes = 20
    bs_days = (binsize_minutes / (60*24))
    min_cycle, max_cycle = np.nanmin(ns.cyclenum_sw), np.nanmax(ns.cyclenum_sw)
    for ix, cycle in enumerate(
        np.arange(min_cycle, max_cycle)
    ):
        sel = (ns.cyclenum_sw == cycle)
        _bd = phase_bin_magseries(ns.x_sw[sel], ns.y_resid_sw[sel],
                                  binsize=bs_days, minbinelems=2)
        color = cmap(norm(cycle))

        if _bd is not None:
            ax.plot(_bd['binnedphases'], gfn(n(_bd['binnedmags'])), zorder=ix,
                    lw=0.3, c=color, alpha=1, rasterized=True)

    ax.plot(ns.x_w_model, n(ns.y_w_model-ns.y_w_model), lw=0.5, zorder=2, c='C0')
    ylim_resid = factor*np.array(get_ylimguess(n(ns.y_resid_sw_b)))
    ax.set_ylim(ylim_resid)
    ax.update({'xlabel': 'φ', 'ylabel': 'Resid [%]', 'xlim':[-0.6, 0.6]})

    # flux vs time raw
    ax = axd['C']
    bd = time_bin_magseries(ns.time, n(ns.flux), binsize=900, minbinelems=1)
    ax.scatter(ns.time, n(ns.flux), s=0.1, marker='o', c='lightgray', linewidths=0)
    ax.scatter(bd['binnedtimes'], bd['binnedmags'], c='k', s=0.3, zorder=2, linewidths=0)
    ylim = get_ylimguess(bd['binnedmags'])
    #ylim = get_ylimguess(n(ns.flux))
    ax.set_ylim(ylim)
    ax.update({'xlabel': '', 'ylabel': '$\Delta$ Flux [%]'})
    ax.set_xticklabels([])

    # flux vs time resid
    ax = axd['D']
    bd = time_bin_magseries(ns.time, n(ns.y_resid_nsnw), binsize=900, minbinelems=1)
    ax.scatter(ns.time, n(ns.y_resid_nsnw), s=0.1, marker='o', c='lightgray', linewidths=0)
    ax.scatter(bd['binnedtimes'], bd['binnedmags'], c='k', s=0.3, zorder=2, linewidths=0)
    ax.set_ylim(ylim)
    ax.update({'xlabel': 'TESS Julian Date [days]', 'ylabel': 'Resid [%]'})

    tform = blended_transform_factory(ax.transData, ax.transAxes)
    tmid = np.nanmean(ns.time)
    t_lo = tmid - d['period']/2
    t_hi = tmid + d['period']/2
    ax.hlines(0.9, t_lo, t_hi, colors='red', alpha=1,
              linestyles='-', zorder=5, linewidths=2, transform=tform)

    fig.tight_layout(h_pad=0)

    # colorbar insanity
    ax = axd['E']
    _p = ax.scatter(_bd['binnedphases'], n(_bd['binnedmags'])+99, zorder=3,
                    s=1, c=cycle*np.ones(len(_bd['binnedmags'])),
                    linewidths=0, cmap=cmap, alpha=1,
                    vmin=min_cycle, vmax=max_cycle, rasterized=True)

    x0,y0,dx,dy = 0.03, 0.2, 0.3, 0.05
    axins1 = inset_axes(ax, width="100%", height="100%",
                        bbox_to_anchor=(x0,y0,dx,dy),
                        loc='lower left',
                        bbox_transform=ax.transAxes)


    cb = fig.colorbar(_p, cax=axins1, orientation="horizontal",
                      extend="neither")
    cb.set_ticks([min_cycle, max_cycle])
    cb.set_ticklabels([int(min_cycle), int(max_cycle)])
    cb.ax.tick_params(labelsize='x-small')
    cb.ax.tick_params(size=0, which='both') # remove the ticks
    cb.ax.yaxis.set_ticks_position('left')
    cb.ax.yaxis.set_label_position('left')

    # initial periodogram
    ax = axd['F']
    d = deepcopy(pgdict)

    ax.plot(d['lsp']['periods'], d['lsp']['lspvals'], c='k', lw=1)
    ax.scatter(d['lsp']['nbestperiods'][:7],
               d['lsp']['nbestlspvals'][:7],
               marker='v', s=5, linewidths=0, edgecolors='none',
               color='k', alpha=0.5, zorder=1000)
    ymin, ymax = ax.get_ylim()
    ax.vlines(d['period'], ymin, ymax, colors='darkgray', alpha=1,
              linestyles='-', zorder=-2, linewidths=1)
    P_harmonics = []
    for ix in range(1,11):
        P_harmonics.append(ix*d['period'])
        P_harmonics.append(d['period']/ix)

    sel = (
        (P_harmonics > np.nanmin(d['lsp']['periods']))
        &
        (P_harmonics < np.nanmax(d['lsp']['periods']))
    )
    ax.vlines(nparr(P_harmonics)[sel], ymin, ymax, colors='darkgray',
              alpha=0.5, linestyles=':', zorder=-2, linewidths=0.5)
    ax.set_ylim([ymin, ymax])
    ax.update({'xlabel': 'Period [d]', 'ylabel': 'PDM Θ', 'xscale': 'log'})

    # initial periodogram text
    ax = axd['I']
    ax.set_axis_off()
    txt = 'Raw periods, $P$\n'
    raw_period = d['lsp']['nbestperiods'][0]
    for _P in d['lsp']['nbestperiods'][:7]:
        txt += f"{_P:.4f} d = {_P*24:.3f} h\n"
    txt_x = 0.
    txt_y = 0.5
    ax.text(txt_x, txt_y, txt, ha='left', va='center', fontsize='x-small', zorder=2,
            transform=ax.transAxes)
    del d

    # residual periodogram
    ax = axd['G']
    d = deepcopy(resid_pgdict)

    ax.plot(d['lsp']['periods'], d['lsp']['lspvals'], c='k', lw=1)
    ax.scatter(d['lsp']['nbestperiods'][:7],
               d['lsp']['nbestlspvals'][:7],
               marker='v', s=5, linewidths=0, edgecolors='none',
               color='k', alpha=0.5, zorder=1000)
    ymin, ymax = ax.get_ylim()
    ax.vlines(d['period'], ymin, ymax, colors='darkgray', alpha=1,
              linestyles='-', zorder=-2, linewidths=1)
    P_harmonics = []
    for ix in range(1,11):
        P_harmonics.append(ix*d['period'])
        P_harmonics.append(d['period']/ix)

    sel = (
        (P_harmonics > np.nanmin(d['lsp']['periods']))
        &
        (P_harmonics < np.nanmax(d['lsp']['periods']))
    )
    ax.vlines(nparr(P_harmonics)[sel], ymin, ymax, colors='darkgray',
              alpha=0.5, linestyles=':', zorder=-2, linewidths=0.5)
    ax.set_ylim([ymin, ymax])
    ax.update({'xlabel': 'Period [d]', 'ylabel': 'Resid PDM Θ', 'xscale': 'log'})

    # residual periodogram text
    ax = axd['J']
    ax.set_axis_off()
    txt = 'Resid periods:\n'
    for _P in d['lsp']['nbestperiods'][:7]:
        txt += f"{_P:.4f} d = {_P*24:.3f} h = {_P/raw_period:.3f} "+"$P$\n"
    txt_x = 0.
    txt_y = 0.5
    ax.text(txt_x, txt_y, txt, ha='left', va='center', fontsize='x-small', zorder=2,
            transform=ax.transAxes)

    #
    # residual phase at next period
    #
    x, y = ns.time, ns.y_resid_nsnw
    period = resid_pgdict['period']
    period2 = resid_pgdict['lsp']['nbestperiods'][1]
    t0 = resid_pgdict['t0']

    binsize_minutes = 10
    bs_days = (binsize_minutes / (60*24))

    _pd_sw = phase_magseries(x, y, period, t0, wrap=1, sort=True)
    r_x_sw = _pd_sw['phase']
    r_y_sw = _pd_sw['mags']
    bd = phase_bin_magseries(r_x_sw, r_y_sw, binsize=bs_days, minbinelems=3)
    r_x_sw_b = bd['binnedphases']
    r_y_sw_b = bd['binnedmags']

    _pd_sw = phase_magseries(x, y, period2, t0, wrap=1, sort=True)
    r_x_sw_2 = _pd_sw['phase']
    r_y_sw_2 = _pd_sw['mags']
    bd = phase_bin_magseries(r_x_sw_2, r_y_sw_2, binsize=bs_days, minbinelems=3)
    r_x_sw_b_2 = bd['binnedphases']
    r_y_sw_b_2 = bd['binnedmags']


    # resid flux vs phase (at best *residual* period)
    ax = axd['H']
    ax.scatter(r_x_sw, n(r_y_sw), zorder=1, s=0.2, c='lightgray',
               linewidths=0)
    ax.scatter(r_x_sw_b, n(r_y_sw_b), zorder=3, s=1, c='k',
               linewidths=0)

    factor = 2
    ylim = factor*np.array(get_ylimguess(n(r_y_sw_b)))
    ax.set_ylim(ylim)
    ax.update({#'xlabel': 'φ',
               'ylabel': 'Resid [%]', 'xlim':[-0.6,0.6],
               'title': f"P$_2$ = {resid_pgdict['period']*24:.3f} hr (best)"
              })

    # resid flux vs phase (at best *residual* period)
    ax = axd['K']

    ax.scatter(r_x_sw_2, n(r_y_sw_2), zorder=1, s=0.2, c='lightgray',
               linewidths=0)
    ax.scatter(r_x_sw_b_2, n(r_y_sw_b_2), zorder=3, s=1, c='k',
               linewidths=0)

    factor = 2
    ylim = factor*np.array(get_ylimguess(n(r_y_sw_b_2)))
    ax.set_ylim(ylim)
    ax.update({#'xlabel': 'φ',
               'ylabel': 'Resid [%]', 'xlim':[-0.6,0.6],
               'title': f"P$_2$ = {period2*24:.3f} hr (2nd)"
              })

    savefig(fig, pngpath, dpi=400, writepdf=0)


def plot_lc_mosaic(outdir, subset_id=None, showtitles=0,
                   titlefontsize='xx-small'):
    """
    this plotter only works on phtess3 or analogous
    """
    # get ticids and sector
    if subset_id in ['dlt150_good_0', 'dlt150_good_1', 'dlt150_good_all']:
        csvpath = join(DATADIR, 'targetlists',
                       '20230411_good_CPV_ticids_d_lt_150pc_sectorpref.csv')
    elif subset_id in ['dlt150_good_changers']:
        csvpath = join(DATADIR, 'targetlists',
                       '20230411_good_CPV_ticids_d_lt_150pc_sectorpref_SIXCHANGERS.csv')
    elif subset_id in ['dlt150_good_allchangers_2count',
                       'dlt150_good_allchangers_3count']:
        csvpath = join(DATADIR, 'targetlists',
                       'before_after_stars_S55_max_2mindataonly_goodCPVonly'
                       '_20230530_20230411_goodandmaybe_CPV_ticids_d_lt_150pc.csv')
    elif subset_id in ['fav3']:
        csvpath = join(DATADIR, 'targetlists', 'fav3.csv')
    else:
        raise NotImplementedError

    df = pd.read_csv(csvpath, comment='#')
    if 'comment' not in df:
        df['comment'] = ''
    sel = ~(df.comment.str.contains("OMIT") == True)
    df = df[sel]

    if subset_id == 'dlt150_good_0':
        df = df[:25]
    elif subset_id == 'dlt150_good_1':
        df = df[25:]
    elif subset_id == 'dlt150_good_changers':
        assert len(df) == 12
    elif subset_id in ['dlt150_good_all', 'fav3']:
        pass
    elif subset_id in ['dlt150_good_allchangers_2count',
                       'dlt150_good_allchangers_3count']:
        from collections import Counter
        r = Counter(df.ticid)
        count_df = pd.DataFrame({'ticid':r.keys(), 'count':r.values()})
        df = df.merge(count_df, how='left', on='ticid')
        if '2count' in subset_id:
            #44 panels
            df = df[df['count'] == 2]
        if '3count' in subset_id:
            #27 panels
            df = df[df['count'] == 3]
    else:
        raise NotImplementedError

    # prepare to get lc data
    LOCAL_DEBUG = 0
    from complexrotators.getters import (
        _get_lcpaths_given_ticid, _get_local_lcpaths_given_ticid
    )
    from complexrotators.lcprocessing import (
        cpv_periodsearch, count_phased_local_minima, prepare_cpv_light_curve
    )
    from cdips.utils import bash_grep

    csvpath = join(SPOCDIR, "gaia_X_spoc2min_merge.csv")
    #gdf = pd.read_csv(csvpath)

    # instantiate plot
    plt.close('all')
    set_style('science')

    factor = 0.8
    if subset_id in ['dlt150_good_0', 'dlt150_good_1']:
        fig, axs = plt.subplots(nrows=5, ncols=5, figsize=(factor*6,factor*7),
                                sharex=True, constrained_layout=True)
    elif subset_id in ['dlt150_good_changers']:
        fig, axs = plt.subplots(nrows=2, ncols=6, figsize=(factor*6,(2.5/7)*factor*7),
                                sharex=True, constrained_layout=True)
    elif subset_id in ['dlt150_good_all']:
        fig, axs = plt.subplots(nrows=5, ncols=8, figsize=(factor*6,(4.5/7)*factor*7),
                                sharex=True, constrained_layout=True)
    elif subset_id in ['fav3']:
        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(1.2*3.5, 1.2*1.25),
                                constrained_layout=True)
    elif subset_id == 'dlt150_good_allchangers_2count':
        fig, axs = plt.subplots(nrows=7, ncols=6, figsize=(7,6),
                                sharex=True, constrained_layout=True)
    elif subset_id == 'dlt150_good_allchangers_3count':
        fig, axs = plt.subplots(nrows=5, ncols=6, figsize=(5,6),
                                sharex=True, constrained_layout=True)
    axs = axs.flatten()

    ix = 0
    for ticid, sector, ax in zip(df.ticid, df.sector, axs):

        #lcdir = f'/nfs/phtess2/ar0/TESS/SPOCLC/sector-{sector}'
        lcdir = '/Users/luke/local/complexrotators/fav3/'
        lcpaths = glob(join(lcdir, f"*{ticid}*_lc.fits"))
        if not len(lcpaths) == 1:
            print('bad lcpaths')
            print(ticid, sector)
            print(lcpaths)
            import IPython; IPython.embed()
            assert 0
        lcpath = lcpaths[0]

        cachedir = '/Users/luke/local/complexrotators/cpv_finding/2023catalog_LGB_RJ_concat'

        # get the relevant light curve data
        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid) = prepare_cpv_light_curve(lcpath, cachedir)

        # get period, t0, and periodogram (PDM or LombScargle)
        d = cpv_periodsearch(
            x_obs, y_flat, starid, cachedir, t0='binmin', periodogram_method='pdm'
        )

        if str(ticid) in ['177309964', '335598085', '234295610']:
            if str(ticid) == '335598085':
                d['period'] *= 0.5
            if str(ticid) == '177309964' and int(sector) == 11:
                d['period'] *= 0.5
            if str(ticid) == '234295610' and int(sector) == 28:
                d['period'] *= 0.5

        bd = time_bin_magseries(d['times'], d['fluxs'], binsize=1200, minbinelems=1)
        ylim = get_ylimguess(1e2*(bd['binnedmags']-np.nanmean(bd['binnedmags'])))

        if showtitles:
            titlestr = f'{ticid}, s{sector}, {d["period"]*24:.1f}h'
        else:
            titlestr = None

        binsize_phase = 1/300

        if subset_id in ['dlt150_good_all']:
            BINMS=1.0
            alpha0=0.15
        else:
            BINMS=1.5
            alpha0=0.3

        plot_phased_light_curve(
            d['times'], d['fluxs'], d['t0'], d['period'], None, ylim=ylim,
            xlim=[-0.6,0.6], binsize_phase=binsize_phase, BINMS=BINMS, titlestr=titlestr,
            showtext=None, showtitle=False, figsize=None, c0='darkgray',
            alpha0=alpha0, c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
            fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
            showxticklabels=False, titlefontsize=titlefontsize
        )

        if not showtitles:
            txt = f'{d["period"]*24:.1f}h'
            tform = ax.transAxes
            props = dict(boxstyle='square', facecolor='white', alpha=0.7, pad=0.15,
                         linewidth=0)
            if subset_id in ['dlt150_good_0', 'dlt150_good_1', 'fav3']:
                fontsize = 'small'
            elif subset_id == 'dlt150_good_all':
                fontsize = 'xx-small'
            else:
                fontsize = 'x-small'

            if subset_id == 'dlt150_good_changers':
                if ix >= 6:
                    ax.text(
                        0.97, 0.05, txt, transform=tform, ha='right',
                        va='bottom', color='k', fontsize=fontsize, bbox=props
                    )
            else:
                ax.text(
                    0.97, 0.05, txt, transform=tform, ha='right',
                    va='bottom', color='k', fontsize=fontsize, bbox=props
                )

        ax.set_xticks([-0.5,0,0.5])

        ylow, yhigh = int(np.ceil(ylim[0]))+1, int(np.floor(ylim[1]))-1
        if np.diff([np.abs(ylow), yhigh]) <= 2:
            ylowabs = np.abs(ylow)
            yhighabs = np.abs(yhigh)
            ylow = -np.min([ylowabs, yhighabs])
            yhigh = np.min([ylowabs, yhighabs])
            if ylow == yhigh == 0:
                # 368129164 led to this
                ylow = -1
                yhigh = 1
            if yhigh >= 10:
                ylow = -9
                yhigh = 9

        ax.set_yticks([ylow, 0, yhigh])
        ax.set_yticklabels([ylow, 0, yhigh])

        id_yticks = {
            "405754448": [-1,0,1],
            "300651846": [-3,0,3],
            "167664935": [-9,0,9],
            "353730181": [-6,0,6],
            '335598085': [-3,0,3],
            '302160226': [-6,0,6],
            '234295610': [-4,0,4],
            '440725886': [-9,0,9],
            '5714469': [-9,0,9],
            '442571495': [-3,0,3],
            '402980664': [-4,0,2],
            'TIC_402980664': [-4,0],
            '141146667': [-7,0,7]
        }
        if str(ticid) in list(id_yticks.keys()):
            yticks = id_yticks[str(ticid)]
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks)
        if str(ticid) == '300651846':
            ax.set_ylim([-5,4])
        if str(ticid) == '353730181':
            ax.set_ylim([-7.5,7.5])
        #if str(ticid) == '402980664':
        #    ax.set_ylim([-2.8, 2.5])
        if str(ticid) == '335598085':
            ax.set_ylim([-4.2, 4.2])
        if str(ticid) == '302160226':
            ax.set_ylim([-9,9])
        if str(ticid) == '167664935':
            ax.set_ylim([-17,17])
        if str(ticid) == '234295610':
            ax.set_ylim([-5.5,4.5])
        if str(ticid) == '440725886':
            ax.set_ylim([-12.5,12.5])
        if str(ticid) == '5714469':
            ax.set_ylim([-11,11])
        if str(ticid) == '442571495':
            ax.set_ylim([-5.5,5.5])
        if str(ticid) == '141146667':
            ax.set_ylim([-11, 11])

        if subset_id == 'dlt150_good_changers':
            if str(ticid) == '404144841':
                ax.set_ylim([-8,7.5])
                ylow, yhigh = -6, 6
            if str(ticid) == '201789285':
                ax.set_ylim([-9.2,9.2])
                ylow, yhigh = -7, 7
            if str(ticid) == '206544316':
                ax.set_ylim([-8,7])
                ylow, yhigh = -6, 6
            if str(ticid) == '234295610':
                ax.set_ylim([-5.2,4.2])
                ylow, yhigh = -4, 4
            if str(ticid) == '224283342':
                ax.set_ylim([-5.5,4])
                ylow, yhigh = -3, 3
            if str(ticid) == '177309964':
                ax.set_ylim([-10,10])
                ylow, yhigh = -6, 6
            if str(ticid) == '2234692':
                ax.set_ylim([-7,7])
                ylow, yhigh = -5, 5
            if str(ticid) == '146539195':
                ax.set_ylim([-6.5,5])
                ylow, yhigh = -4, 4
            ax.set_yticks([ylow, 0, yhigh])
            ax.set_yticklabels([ylow, 0, yhigh])

        if subset_id in ['dlt150_good_0', 'dlt150_good_1', 'fav3']:
            labelsize = 'small'
        elif subset_id == 'dlt150_good_all':
            labelsize = 'xx-small'
        else:
            labelsize = 'x-small'

        if subset_id == 'dlt150_good_all':
            ax.minorticks_off()

        if subset_id == 'dlt150_good_all':
            ax.tick_params(axis='both', which='major', labelsize=labelsize,
                           pad=1.5)
        else:
            ax.tick_params(axis='both', which='major', labelsize=labelsize)

        ix += 1

    for ax in axs[-4:]:
        ax.set_xticklabels(['-0.5','0','0.5'])

    if subset_id == 'dlt150_good_1':
        axs[-1].set_axis_off()

    if subset_id in ['dlt150_good_0', 'dlt150_good_1']:
        fs = 'large'
    else:
        fs = 'medium'

    if subset_id != 'fav3':
        fig.text(0.5,-0.01, r"Phase, φ", fontsize=fs)
        fig.text(-0.02,0.5, "$\Delta$ Flux [%]", va='center', rotation=90, fontsize=fs)
    else:
        axs[0].set_ylabel("$\Delta$ Flux [%]", fontsize=fs)
        axs[1].set_xlabel(r"Phase, φ", fontsize=fs)

    # set naming options
    s = f'{subset_id}'
    if showtitles:
        s += '_showtitles'

    fig.tight_layout(h_pad=0.2, w_pad=0.)

    outpath = join(outdir, f'lc_mosaic_{s}.png')
    savefig(fig, outpath, dpi=400)


def plot_full_lcmosaic(outdir, showtitles=1, titlefontsize=3.5,
                       sortby=None, binarytitle=0, rasterized=False):
    """
    the one that makes the publication
    (on mico, assuming that find_CPVs.py has been run using sample_id=="2023catalog_LGB_RJ_concat")
    """

    # load ticids
    N_objects = 50 # fine if longer/shorter, just requires tweaking

    tablepath = join(
        TABLEDIR, "2023_catalog_table",
        '20230613_LGB_RJ_CPV_TABLE_supplemental_selfnapplied.csv'
    )
    df = pd.read_csv(tablepath, sep="|")
    _df = df[df.quality == 1]
    assert len(_df) == N_objects

    # get manually written sectors to plot
    sectorpath = join(
        DATADIR, 'targetlists', '20230613_LGB_RJ_sectorpref.csv'
    )
    tdf = pd.read_csv(sectorpath)

    # merge to same dataframe
    df = _df.merge(tdf, on='ticid', how='inner')
    assert len(df) == N_objects

    if isinstance(sortby, str):
        if sortby in ['tic8_Tmag','tlc_mean_period','dist_pc']:
            df = df.sort_values(by=sortby)
        elif sortby in ['Rcr_over_Rstar']:
            df = df.sort_values(by=sortby, ascending=False)
        elif sortby == 'Ngoodsectors_tic8_Tmag':
            fn = lambda x: len(x.split(','))
            df['Ngoodsectors'] = df.goodsectors_x.apply(fn)
            df = df.sort_values(
                by=['Ngoodsectors', 'tic8_Tmag'],
                ascending=[False, True]
            )
        else:
            raise NotImplementedError

    if binarytitle:
        csvpath = join(TABLEDIR, "2023_catalog_table",
                       "20230613_LGB_RJ_CPV_TABLE_supplemental_selfnapplied.csv")
        tabdf = pd.read_csv(csvpath, sep="|", dtype={'binarityflag':str})

    # prepare to get lc data
    from complexrotators.lcprocessing import (
        cpv_periodsearch, prepare_cpv_light_curve
    )

    #
    # instantiate plot
    #
    plt.close('all')
    set_style('science')

    factor = 0.8
    ncols = 5
    fig, axs = plt.subplots(nrows=10, ncols=ncols,
                            figsize=(factor*7,factor*8.7),
                            constrained_layout=True)
    axs = axs.flatten()

    ix = 0

    for ticid, sector, ax in zip(df.ticid, df.showsector, axs):

        lcdir = '/Users/luke/.lightkurve/cache/mastDownload/TESS'
        lcpaths = glob(join(
            lcdir, f"tess*s{str(sector).zfill(4)}*{ticid}*",
            f"*{ticid}*_lc.fits")
        )
        if not len(lcpaths) == 1:
            print('bad lcpaths')
            print(ticid, sector)
            print(lcpaths)
            import IPython; IPython.embed()
            assert 0

        lcpath = lcpaths[0]

        cachedir = '/Users/luke/local/complexrotators/cpv_finding/2023catalog_LGB_RJ_concat_BACKUP'

        # get the relevant light curve data
        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid) = prepare_cpv_light_curve(lcpath, cachedir)

        # get period, t0, and periodogram (PDM or LombScargle)
        d = cpv_periodsearch(
            x_obs, y_flat, starid, cachedir, t0='binmin', periodogram_method='pdm'
        )

        if str(ticid) in ['177309964', '335598085', '234295610']:
            if str(ticid) == '335598085':
                d['period'] *= 0.5
            if str(ticid) == '177309964' and int(sector) == 11:
                d['period'] *= 0.5
            if str(ticid) == '234295610' and int(sector) == 28:
                d['period'] *= 0.5

        bd = time_bin_magseries(d['times'], d['fluxs'], binsize=1200, minbinelems=1)

        ylim = get_ylimguess(1e2*(bd['binnedmags']-np.nanmean(bd['binnedmags'])))

        if showtitles:
            titlestr = f'TIC{ticid}, S{sector}, {d["period"]*24:.1f}h'
            if binarytitle:
                binflag = tabdf.loc[
                    tabdf.ticid.astype(str)==str(ticid), 'binarityflag'
                ]
                titlestr += f", {binflag.iloc[0]}"
        else:
            titlestr = None

        binsize_phase = 1/300

        BINMS=1.0
        alpha0=0.15

        plot_phased_light_curve(
            d['times'], d['fluxs'], d['t0'], d['period'], None, ylim=ylim,
            xlim=[-0.6,0.6], binsize_phase=binsize_phase, BINMS=BINMS, titlestr=titlestr,
            showtext=None, showtitle=False, figsize=None, c0='darkgray',
            alpha0=alpha0, c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
            fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
            showxticklabels=False, titlefontsize=titlefontsize,
            titlepad=0.05, rasterized=rasterized
        )

        if not showtitles:
            txt = f'{d["period"]*24:.1f}h'
            tform = ax.transAxes
            props = dict(boxstyle='square', facecolor='white', alpha=0.7, pad=0.15,
                         linewidth=0)
            fontsize = 'xx-small'
            ax.text(
                0.97, 0.05, txt, transform=tform, ha='right',
                va='bottom', color='k', fontsize=fontsize, bbox=props
            )

        ax.set_xticks([-0.5,0,0.5])

        ylow, yhigh = int(np.ceil(ylim[0]))+1, int(np.floor(ylim[1]))-1
        if np.diff([np.abs(ylow), yhigh]) <= 2:
            ylowabs = np.abs(ylow)
            yhighabs = np.abs(yhigh)
            ylow = -np.min([ylowabs, yhighabs])
            yhigh = np.min([ylowabs, yhighabs])
            if ylow == yhigh == 0:
                # 368129164 led to this
                ylow = -1
                yhigh = 1
            if yhigh >= 10:
                ylow = -9
                yhigh = 9

        ax.set_yticks([ylow, 0, yhigh])
        ax.set_yticklabels([ylow, 0, yhigh])

        id_yticks = {
            "405754448": [-1,0,1],
            "300651846": [-3,0,3],
            "167664935": [-9,0,9],
            "353730181": [-6,0,6],
            '335598085': [-3,0,3],
            '302160226': [-3,0,2],
            '234295610': [-4,0,4],
            '440725886': [-9,0,9],
            '5714469': [-6,0,6],
            '442571495': [-3,0,3],
            '402980664': [-2,0,1],
            '141146667': [-7,0,7],
            '363963079': [-7,0,7],
            '89463560': [-4,0,4],
            '142173958': [-8,0,8],
            '405910546': [-6,0,3],
            '118449916': [-6,0,4],
            '425937691': [-8,0,8],
            '397791443': [-9,0,9],
            '312410638': [-4,0,4],
            '332517282': [-7,0,7],
            '38539720': [-8,0,8],
        }
        if str(ticid) in list(id_yticks.keys()):
            yticks = id_yticks[str(ticid)]
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks)
            if yticks[2] >= 9:
                diff = 4
            elif yticks[2] > 5:
                diff = 2.5
            elif yticks[2] > 2:
                diff = 1.5
            else:
                diff = 1.
            ax.set_ylim([yticks[0]-diff, yticks[2]+diff])

            if str(ticid) == '397791443':
                ax.set_ylim([-17,13])

        labelsize = 'xx-small'

        ax.minorticks_off()

        ax.tick_params(axis='both', which='major', labelsize=labelsize,
                       pad=1.5)

        ix += 1

    for ax in axs[-ncols:]:
        ax.set_xticklabels(['-0.5','0','0.5'])
    #axs[-1].set_axis_off()

    fs = 'medium'
    fig.text(0.5,0.0, r"Phase, φ", fontsize=fs)
    fig.text(-0.01,0.5, "$\Delta$ Flux [%]", va='center', rotation=90, fontsize=fs)

    # set naming options
    s = ''
    if showtitles:
        s += '_showtitles'
    if binarytitle:
        s += '_binarytitle'
    if isinstance(sortby, str):
        s += f'_{sortby}'
    if rasterized:
        s += '_rasterized'

    # height/width
    fig.tight_layout(h_pad=0.06, w_pad=0.)

    #axs[-1].set_axis_off()

    outpath = join(outdir, f'full_lcmosaic{s}.png')
    savefig(fig, outpath, dpi=600)


def plot_beforeafter_mosaic(outdir, showtitles=1, titlefontsize=3.75,
                            sortby=None, rasterized=False, fav3=0):
    """
    the one that makes the publication
    (on mico, assuming that find_CPVs.py has been run using
    sample_id=="2023catalog_LGB_RJ_concat")
    """

    # load ticids
    N_objects = 50 # fine if longer/shorter, just requires tweaking
    if fav3:
        N_objects = 3

    tablepath = join(
        TABLEDIR, "2023_catalog_table",
        '20230613_LGB_RJ_CPV_TABLE_supplemental_selfnapplied.csv'
    )
    df = pd.read_csv(tablepath, sep="|")
    _df = df[~pd.isnull(df.goodsectors)]

    if fav3:
        ticids = np.array(['402980664', '206544316', '404144841'])
        _df = _df[_df.ticid.astype(str).isin(ticids)]

    assert len(_df) == N_objects

    # get manually written sectors to plot
    sectorpath = join(
        DATADIR, 'targetlists', '20230613_LGB_RJ_changerspref.csv'
    )
    tdf = pd.read_csv(sectorpath, comment='#')

    # merge to same dataframe
    df = _df.merge(tdf, on='ticid', how='inner')
    assert len(df) == N_objects

    if isinstance(sortby, str):
        if sortby in ['tic8_Tmag','tlc_mean_period','dist_pc']:
            df = df.sort_values(by=sortby)
        elif sortby == 'Ngoodsectors_tic8_Tmag':
            fn = lambda x: len(x.split(','))
            df['Ngoodsectors'] = df.goodsectors_x.apply(fn)
            df = df.sort_values(
                by=['Ngoodsectors', 'tic8_Tmag'],
                ascending=[False, True]
            )
        else:
            raise NotImplementedError

    # require "changersectors" manual label to be filled -- these are all cases
    # for which "sectors" had >26 day baselines (irrespective of 2-day-data
    # availability)
    df = df[~pd.isnull(df.changersectors)]

    fn = lambda x: x.split('-')[0]
    df['beforesector'] = df['changersectors'].apply(fn)
    fn = lambda x: x.split('-')[1]
    df['aftersector'] = df['changersectors'].apply(fn)

    #
    # attempt pre-download all 120-second cadence light curves.  both the
    # before and after sectors have to exist at 120-second cadence.
    #
    from complexrotators.getters import _get_lcpaths_fromlightkurve_given_ticid
    beforesectorfound, aftersectorfound = [], []
    for sectorcol in ['beforesector','aftersector']:
        for ticid, sector in zip(df['ticid'], df[sectorcol]):
            lcdir = '/Users/luke/.lightkurve/cache/mastDownload/TESS'
            lcpaths = glob(join(
                lcdir, f"tess*s{str(sector).zfill(4)}*{ticid}*",
                f"*{ticid}*_lc.fits")
            )
            if not len(lcpaths) == 1:
                lcpaths = _get_lcpaths_fromlightkurve_given_ticid(str(ticid), require_lc=1)
                thepath = [l for l in lcpaths if f's{str(sector).zfill(4)}' in l]
                if len(thepath) == 1:
                    if sectorcol == 'beforesector':
                        beforesectorfound.append(True)
                    elif sectorcol == 'aftersector':
                        aftersectorfound.append(True)
                else:
                    if sectorcol == 'beforesector':
                        beforesectorfound.append(False)
                    elif sectorcol == 'aftersector':
                        aftersectorfound.append(False)
            else:
                if sectorcol == 'beforesector':
                    beforesectorfound.append(True)
                elif sectorcol == 'aftersector':
                    aftersectorfound.append(True)

    df['beforesectorfound'] = np.array(beforesectorfound)
    df['aftersectorfound'] = np.array(aftersectorfound)

    sel = df['beforesectorfound'] & df['aftersectorfound']

    df = df[sel]

    N = len(df)
    print(f"Got {N} stars with beforesectorfound and aftersectorfound...")
    cut_N = 27 # 54 panels total
    if fav3:
        cut_N = 3
    print(f"Cutting to 27...")
    df = df.head(n=cut_N)

    a_ticids, a_beforesectors, a_aftersectors = (
        np.array(df.ticid), np.array(df.beforesector), np.array(df.aftersector)
    )
    ticids, sectors = [], []
    for ix in range(cut_N*2):
        ticids.append(a_ticids[ix // 2])
        if ix % 2 == 0:
            sectors.append(a_beforesectors[ix // 2])
        else:
            sectors.append(a_aftersectors[ix // 2])

    # prepare to get lc data
    from complexrotators.lcprocessing import (
        cpv_periodsearch, prepare_cpv_light_curve
    )

    #
    # instantiate plot
    #
    plt.close('all')
    set_style('science')


    SIMPLEGRID = 0

    if SIMPLEGRID or fav3:
        if not fav3:
            factor = 0.8
            # this gets it looking 90% correct
            fig, axs = plt.subplots(nrows=9, ncols=6,
                                    figsize=(factor*7,factor*8.7),
                                    constrained_layout=True)
            axs = axs.flatten()
        else:
            factor = 0.8
            # this gets it looking 90% correct
            fig, axs = plt.subplots(nrows=3, ncols=2,
                                    figsize=(factor*2,factor*3),
                                    constrained_layout=True)
            axs = axs.flatten()

    else:
        # this is the correct way to do it
        # see tests/gridspec_demo3.py, where this insanity was prototyped
        factor = 0.75
        fig = plt.figure(figsize=(factor*6.5,factor*9.5))

        from matplotlib.gridspec import GridSpec

        x0 = 0.03
        smx = 0.03
        dx = 0.30
        hspace = 0.2
        wspace = 0.02
        left, right = x0, x0+dx

        axs = []

        gs1 = GridSpec(9, 2, left=left, right=right, wspace=wspace, hspace=hspace)
        for ix in range(9):
            ax = fig.add_subplot(gs1[ix, 0])
            axs.append(ax)
            ax = fig.add_subplot(gs1[ix, 1])
            axs.append(ax)

        left, right = x0+smx+dx, x0+smx+2*dx
        gs2 = GridSpec(9, 2, left=left, right=right, wspace=wspace, hspace=hspace)
        for ix in range(9):
            ax = fig.add_subplot(gs2[ix, 0])
            axs.append(ax)
            ax = fig.add_subplot(gs2[ix, 1])
            axs.append(ax)

        left, right = x0+2*smx+2*dx, x0+2*smx+3*dx
        gs3 = GridSpec(9, 2, left=left, right=right, wspace=wspace, hspace=hspace)
        for ix in range(9):
            ax = fig.add_subplot(gs3[ix, 0])
            axs.append(ax)
            ax = fig.add_subplot(gs3[ix, 1])
            axs.append(ax)


    ix = 0

    _periods = []
    for ticid, sector, ax in zip(ticids, sectors, axs):

        lcdir = '/Users/luke/.lightkurve/cache/mastDownload/TESS'
        lcpaths = glob(join(
            lcdir, f"tess*s{str(sector).zfill(4)}*{ticid}*",
            f"*{ticid}*_lc.fits")
        )
        if not len(lcpaths) == 1:
            print('bad lcpaths')
            print(ticid, sector)
            print(lcpaths)
            import IPython; IPython.embed()
            assert 0

        lcpath = lcpaths[0]

        cachedir = '/Users/luke/local/complexrotators/cpv_finding/2023catalog_LGB_RJ_concat'

        # get the relevant light curve data
        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid) = prepare_cpv_light_curve(lcpath, cachedir)

        # get period, t0, and periodogram (PDM or LombScargle)
        d = cpv_periodsearch(
            x_obs, y_flat, starid, cachedir, t0='binmin', periodogram_method='pdm'
        )

        if str(ticid) == '335598085':
            d['period'] = 1.3211004981093784/2
        if str(ticid) == '177309964' and int(sector) == 11:
            d['period'] *= 0.5
        if str(ticid) == '234295610' and int(sector) == 28:
            d['period'] *= 0.5
        if str(ticid) == '272248916':
            d['period'] = 0.3708350878122453
        if str(ticid) == '289840926':
            d['period'] = 0.1999953362691368 # lol what

        _periods.append(d['period'])

        bd = time_bin_magseries(d['times'], d['fluxs'], binsize=1200, minbinelems=1)

        ylim = get_ylimguess(1e2*(bd['binnedmags']-np.nanmean(bd['binnedmags'])))

        if showtitles:
            if SIMPLEGRID:
                if ix % 2 == 0:
                    titlestr = f'TIC{ticid}, S{sector}, {d["period"]*24:.1f}h'
                else:
                    titlestr = f'... S{sector}'
            else:
                if ix % 2 == 0:
                    spstr = 49*' '
                    titlestr = spstr+f'TIC{ticid}, {d["period"]*24:.1f}h'
                else:
                    titlestr = None
        else:
            titlestr = None

        VERBOSE = 1
        if VERBOSE:
            print(ticid, ix, f'{d["period"]*24:.1f}h')

        binsize_phase = 1/300

        BINMS=1.0
        alpha0=0.15

        # visually distinct colors generated by https://mokole.com/palette.html
        # 27 colors; 5^ minimum lum, 80% max, 5000 loops
        # sub indigo for lime, peru for gold, and darkorange for khaki
        if not fav3:
            _colors = [
                "#696969",
                "#8b4513",
                "#006400",
                "#808000",
                "#483d8b",
                "#3cb371",
                "#4682b4",
                "#000080",
                "#9acd32",
                "#8b008b",
                "#b03060",
                "#48d1cc",
                "#ff4500",
                "#ff8c00",
                #"#ffd700", # gold
                "peru",
                #"#00ff00", # lime
                "indigo",
                "#8a2be2",
                "#00ff7f",
                "#dc143c",
                "#0000ff",
                "#d8bfd8",
                "#ff00ff",
                "#1e90ff",
                #"#f0e68c", # khaki
                "darkorange",
                "#ff1493",
                "#ffa07a",
                "#ee82ee",
            ]
            np.random.seed(42)
            np.random.shuffle(_colors)
            #if ix % 2 == 0:
            #    c1 = f"C{ix % 10}"
            #else:
            #    c1 = f"C{(ix-1) % 10}"
        elif fav3:
            _colors = [f"C{ix}" for ix in range(3)]
        else:
            raise NotImplementedError

        c1 = _colors[ix // 2]

        if fav3:
            titlestrcopy = deepcopy(titlestr)
            titlestr = ''

        plot_phased_light_curve(
            d['times'], d['fluxs'], d['t0'], d['period'], None, ylim=ylim,
            xlim=[-0.6,0.6], binsize_phase=binsize_phase, BINMS=BINMS, titlestr=titlestr,
            showtext=None, showtitle=False, figsize=None, c0='darkgray',
            alpha0=alpha0, c1=c1, alpha1=1, phasewrap=True, plotnotscatter=False,
            fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
            showxticklabels=False, titlefontsize=titlefontsize,
            titlepad=0.05, rasterized=rasterized
        )

        if SIMPLEGRID:
            pass
        else:
            txt = f"S{sector}"
            props = dict(boxstyle='square', facecolor='white', alpha=0.7,
                         pad=0.15, linewidth=0)
            tform = ax.transAxes
            if ix % 2 == 0:
                # lower right
                ax.text(
                    0.93, 0.09, txt, transform=tform, ha='right', va='bottom',
                    color='k', fontsize=titlefontsize, bbox=props, zorder=99
                )
            else:
                # lower left
                ax.text(
                    0.07, 0.09, txt, transform=tform, ha='left', va='bottom',
                    color='k', fontsize=titlefontsize, bbox=props, zorder=99
                )

        if not showtitles:
            txt = f'{d["period"]*24:.1f}h'
            tform = ax.transAxes
            props = dict(boxstyle='square', facecolor='white', alpha=0.7, pad=0.15,
                         linewidth=0)
            fontsize = 'xx-small'
            ax.text(
                0.97, 0.05, txt, transform=tform, ha='right',
                va='bottom', color='k', fontsize=fontsize, bbox=props
            )

        ax.set_xticks([-0.5,0,0.5])

        ylow, yhigh = int(np.ceil(ylim[0]))+1, int(np.floor(ylim[1]))-1
        if np.diff([np.abs(ylow), yhigh]) <= 2:
            ylowabs = np.abs(ylow)
            yhighabs = np.abs(yhigh)
            ylow = -np.min([ylowabs, yhighabs])
            yhigh = np.min([ylowabs, yhighabs])
            if ylow == yhigh == 0:
                # 368129164 led to this
                ylow = -1
                yhigh = 1
            if yhigh >= 10:
                ylow = -9
                yhigh = 9

        ax.set_yticks([ylow, 0, yhigh])
        ax.set_yticklabels([ylow, 0, yhigh], fontsize='small')

        id_yticks = {
            "405754448": [-1,0,1],
            "300651846": [-3,0,3],
            "167664935": [-9,0,9],
            "353730181": [-6,0,6],
            '335598085': [-3,0,3],
            '302160226': [-3,0,2],
            '234295610': [-4,0,4],
            '201898222': [-2,0,2],
            '264767454': [-1,0,1],
            '440725886': [-9,0,9],
            '5714469': [-6,0,6],
            '442571495': [-3,0,3],
            '264599508': [-3,0,3],
            '193831684': [-9,0,4],
            '177309964': [-7,0,7],
            '402980664': [-2,0,1],
            '141146667': [-7,0,7],
            '363963079': [-7,0,7],
            '89463560': [-4,0,4],
            '142173958': [-7,0,7],
            '405910546': [-6,0,3],
            '118449916': [-6,0,4],
            '425937691': [-8,0,8],
            '397791443': [-9,0,9],
            '312410638': [-4,0,4],
            '332517282': [-7,0,7],
            '144486786': [-1,0,1],
            '38820496': [-1,0,1],
            '38539720': [-8,0,8],
            '289840926': [-9,0,9],
            '404144841': [-6,0,6],
            '89463560': [-4,0,4],
            '368129164': [-1,0,1],
            '50745567': [-2,0,2],
            '425933644': [-4,0,4],
            '146539195': [-4,0,3],
            '206544316': [-6,0,6],
            '272248916': [-2,0,2],
            '178155030': [-3,0,3],
            '224283342': [-4,0,2],

        }
        if str(ticid) in list(id_yticks.keys()):
            yticks = id_yticks[str(ticid)]
            ax.set_yticks(yticks)

            if ix % 2 == 0:
                ax.set_yticklabels(yticks)
            else:
                ax.set_yticklabels([])

            if yticks[2] >= 9:
                diff = 4
            elif yticks[2] > 5:
                diff = 2.5
            elif yticks[2] > 2:
                diff = 1.5
            else:
                diff = 1.
            ax.set_ylim([yticks[0]-diff, yticks[2]+diff])

            if str(ticid) == '397791443':
                ax.set_ylim([-17,13])
            if str(ticid) == '402980664':
                ax.set_ylim([-3.5,2.5])
            if str(ticid) == '206544316' and fav3:
                ax.set_ylim([-8,6])
                ax.set_yticks([-6,0,4])
                if ix % 2 == 0:
                    ax.set_yticklabels([-6,0,4])

        labelsize = 'xx-small'

        ax.minorticks_off()

        ax.tick_params(axis='both', which='major', labelsize=labelsize,
                       pad=1.5)

        ix += 1

    if not fav3:
        outdf = pd.DataFrame({
            'ticid': ticids,
            'sector': sectors,
            'period': _periods
        })
        outcsv = join(outdir, 'beforeafter_periods.csv')
        outdf.to_csv(outcsv, index=False)

    # set x tick labels
    if SIMPLEGRID and not fav3:
        for ax in axs[-6:]:
            ax.set_xticklabels(['-0.5','0','0.5'], fontsize='small')
    elif fav3:
        for ix in [4,5]:
            axs[ix].set_xticklabels(['-0.5','0','0.5'], fontsize=labelsize)
    else:
        for ix in [16,17, 34,35, 52,53]:
            axs[ix].set_xticklabels(['-0.5','0','0.5'], fontsize=labelsize)

    fs = 'medium'
    if fav3:
        fs = 'small'
    if SIMPLEGRID:
        fig.text(0.5,0.0, r"Phase, φ", fontsize=fs)
        fig.text(-0.01,0.5, "$\Delta$ Flux [%]", va='center', rotation=90, fontsize=fs)
    elif fav3:
        fig.text(0.5,0.0, r"Phase, φ", fontsize=fs, va='bottom', ha='center')
        fig.text(0.02,0.5, "$\Delta$ Flux [%]", va='center', ha='center', rotation=90, fontsize=fs)
    else:
        fig.text(0.5,0.07, r"Phase, φ", fontsize=fs, va='bottom', ha='center')
        fig.text(-0.01,0.5, "$\Delta$ Flux [%]", va='center', ha='center', rotation=90, fontsize=fs)


    # set naming options
    s = ''
    if showtitles:
        s += '_showtitles'
    if isinstance(sortby, str):
        s += f'_{sortby}'
    if rasterized:
        s += '_rasterized'
    if fav3:
        s += '_fav3'

    # height/width
    if SIMPLEGRID:
        fig.tight_layout(h_pad=0.1, w_pad=0.)
    if fav3:
        # doctor it after
        fig.tight_layout(h_pad=1.1, w_pad=0.2)

    outpath = join(outdir, f'beforeafter_mosaic{s}.png')
    savefig(fig, outpath, dpi=600)



def plot_cadence_comparison(outdir, ticid: str = None, sector: int = None):

    assert isinstance(ticid, str)
    assert isinstance(sector, int)

    # get data
    lc_cadences = '2min'
    lclist = _get_cpv_lclist(lc_cadences, "TIC "+ticid)

    if len(lclist) == 0:
        print(f'WRN! Did not find light curves for {ticid}. Escaping.')
        return 0

    # for each light curve (sector / cadence specific), detrend if needed, get
    # the best period, and then phase-fold.
    d = None

    for lc in lclist:

        if lc.sector != sector:
            continue

        (time, flux, qual, x_obs, y_obs, y_flat,
         y_trend, x_trend, cadence_sec, sector,
         starid) = prepare_given_lightkurve_lc(lc, ticid, outdir)

        # get t0, period, lsp
        d = cpv_periodsearch(x_obs, y_flat, starid, outdir, t0='binmin')

    if d is None:
        print(f'WRN! Missed LC for {ticid}, SECTOR {sector}. Escaping.')
        return 0

    bd_1800 = time_bin_magseries(
        x_obs, y_obs,
        binsize=1800, minbinelems=1
    )
    bd_600 = time_bin_magseries(
        x_obs, y_obs,
        binsize=600, minbinelems=1
    )
    #bd_1800 = phase_bin_magseries(
    #    bd_1800['binnedtimes'], bd_1800['binnedmags'],
    #    binsize=0.005, minbinelems=1
    #)
    #bd_600 = phase_bin_magseries(
    #    bd_600['binnedtimes'], bd_600['binnedmags'],
    #    binsize=0.005, minbinelems=1
    #)
    #bd_120 = phase_bin_magseries(
    #    d['times'], d['fluxs'],
    #    binsize=0.005, minbinelems=1
    #)

    # make plot
    plt.close('all')
    set_style("clean")

    fig, ax = plt.subplots(figsize=(2,2))

    plot_t0 = d['t0']
    plot_period = d['period']

    plot_phased_light_curve(
        bd_1800['binnedtimes'], bd_1800['binnedmags'], plot_t0, plot_period, None,
        fig=fig, ax=ax,
        binsize_phase=0.01, xlim=[-0.6,0.6], yoffset=20, showtext='1800 sec',
        savethefigure=False, dy=-5, xtxtoffset=-0.24
    )
    plot_phased_light_curve(
        bd_600['binnedtimes'], bd_600['binnedmags'], plot_t0, plot_period, None,
        fig=fig, ax=ax,
        binsize_phase=0.01, xlim=[-0.6,0.6], yoffset=10, showtext='600 sec',
        savethefigure=False, dy=-5, xtxtoffset=-0.24
    )
    plot_phased_light_curve(
        x_obs, y_obs, plot_t0, plot_period, None,
        fig=fig, ax=ax,
        binsize_phase=0.01, xlim=[-0.6,0.6], yoffset=0, showtext='120 sec',
        savethefigure=False, dy=-5, xtxtoffset=-0.24
    )
    ax.set_xlabel('Phase, φ')
    ax.set_ylabel('$\Delta$ Flux [%]')
    ax.set_ylim([-9,28])

    ax.set_yticks([-5, 5, 15, 25])
    ax.set_yticklabels([-5, 5, 15, 25])

    ax.set_xticks([-0.5, 0, 0.5])
    ax.set_xticklabels([-0.5, 0, 0.5])

    # set naming options
    s = f'{ticid}_S{sector}'

    outpath = join(outdir, f'cadence_comparison_{s}.png')
    savefig(fig, outpath, dpi=400)


def plot_tic4029_segments(outdir):

    plt.close('all')
    set_style('science')

    f = 0.7
    fig, ax = plt.subplots(figsize=(f*7.5, f*9))

    ############
    # get data #
    ############
    ticid = '402980664'
    sample_id = '2023catalog_LGB_RJ_concat' # used for cacheing
    period = 18.5611/24
    t0 = 1791.12
    cachedir = join(LOCALDIR, "cpv_finding", sample_id)

    from complexrotators.getters import _get_lcpaths_fromlightkurve_given_ticid
    from complexrotators.lcprocessing import prepare_cpv_light_curve

    lcpaths = _get_lcpaths_fromlightkurve_given_ticid(ticid)

    # stack over all sectors
    times, fluxs, cadencenos = [], [], []
    for lcpath in np.sort(lcpaths):
        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid, cadenceno) = prepare_cpv_light_curve(
             lcpath, cachedir, returncadenceno=1
         )
        times.append(x_obs)
        fluxs.append(y_flat)
        cadencenos.append(cadenceno)
    times = np.hstack(np.array(times, dtype=object).flatten())
    fluxs = np.hstack(np.array(fluxs, dtype=object).flatten())
    cadencenos = np.hstack(np.array(cadencenos, dtype=object).flatten())

    # write copy of LC for glue cadence flagging
    csvpath = join(outdir, f"stitched_TIC{ticid}_t_f_cadence.csv")
    if not os.path.exists(csvpath):
        _df = pd.DataFrame({'time':times, 'flux':fluxs, 'cadenceno':cadencenos})
        _df.to_csv(csvpath, index=False)
        print(f"Wrote {csvpath}")
    else:
        print(f"Found {csvpath}")

    from astrobase.lcmath import find_lc_timegroups
    ngroups, groups = find_lc_timegroups(times, mingap=12/24)

    normfunc = True

    for n, g in enumerate(groups):

        gtime = times[g]
        gflux = fluxs[g]

        if normfunc:
            gflux = gflux - np.nanmean(gflux)

        # t = t0 + P*e
        e_start = int(np.floor((gtime[0] - t0)/period))
        e_end = int(np.floor((gtime[-1] - t0)/period))

        if e_end - e_start < 1:
            continue

        if e_start == -1:
            e_start = 0

        txt = f"{e_start}"+"$\,$-$\,$"+f"{e_end}"

        bd = time_bin_magseries(gtime, gflux, binsize=900, minbinelems=1)

        if e_start > 1000:
            # omitted segment hack
            n -= 1

        yoffset = -6 * n

        norm = lambda _y: 1e2*_y + yoffset

        y = norm(bd['binnedmags'])
        x = bd['binnedtimes']
        x0 = np.nanmin(x)
        x -= x0

        sx, sy, _ = sigclip_magseries(x, y, np.zeros_like(y), sigclip=[10,2.5],
                                      iterative=False, niterations=None,
                                      meanormedian='median',
                                      magsarefluxes=True)

        # show flares in light gray; everything else in black; don't draw lines
        # between time gaps.
        _, _groups = find_lc_timegroups(x, mingap=0.5/24)
        for _g in _groups:
            ax.plot(x[_g], y[_g], c='k', zorder=1, lw=0.2, alpha=0.25)

        _, _groups = find_lc_timegroups(sx, mingap=0.5/24)
        for _g in _groups:
            ax.plot(sx[_g], sy[_g], c='k', zorder=2, lw=0.3)

        # manage text
        fontsize = 6
        tform = blended_transform_factory(ax.transAxes, ax.transData)
        props = dict(boxstyle='square', facecolor='white', alpha=0.7, pad=0.15,
                     linewidth=0)
        ax.text(0.97,
                np.nanpercentile(norm(bd['binnedmags']), 98), txt,
                transform=tform, ha='right',va='bottom', color='k',
                fontsize=fontsize, bbox=props)

    x, y, dx, dy = 7.52, -44, 0, 2
    ax.vlines(x, y, y+dy, colors='red', zorder=30, linewidths=1.4)

    x, y, dx, dy = 5.98, -51, 0, 2
    ax.vlines(x, y, y+dy, colors='red', zorder=30, linewidths=1.4)

    x, y, dx, dy = 10.05, -26.5, 0, 2
    ax.vlines(x, y, y+dy, colors='red', zorder=30, linewidths=1.4)


    ax.text(0.97, 0.98, 'Cycle #', va='top', ha='right', fontsize=7,
            bbox=props, transform=ax.transAxes)

    ax.set_xlabel('Days since segment began')
    ax.set_ylabel('$\Delta$ Flux [%]')

    ax.tick_params(axis='both', which='major', labelsize='small')

    fig.tight_layout()

    # set naming options
    s = ''

    bn = 'tic4029_segments'
    outpath = join(outdir, f'{bn}{s}.png')
    savefig(fig, outpath, dpi=500)


def plot_catalogscatter(outdir, showmaybe=0, plotsubset=None, emphruwe=1):

    # get data
    tablepath = join(
        TABLEDIR, "2023_catalog_table",
        '20230613_LGB_RJ_CPV_TABLE_supplemental_selfnapplied.csv'
    )
    df = pd.read_csv(tablepath, sep="|")
    # drop "debunked" CQVs
    maybe_df = df[df.quality == 0]
    df = df[df.quality == 1]
    binary_df = pd.concat((
        df[
            ( (df.dr3_ruwe > 2) & (df.multipleperiodflag.astype(int)==1) )
            |
            (df.rvscatterflag.astype(int)==1)
        ],
        maybe_df[
            ((maybe_df.dr3_ruwe > 2) & (maybe_df.multipleperiodflag.astype(int)==1))
            |
            (maybe_df.rvscatterflag.astype(int)==1)
        ]
    ))

    #y, x, ylabel, xlabel, yscale, xscale
    if plotsubset is None:
        tuples = [
            ("M_G", "bp_rp", '$\mathrm{M}_{G}$ [mag]', '$G_{\mathrm{BP}}-G_{\mathrm{RP}}$ [mag]', 'linear', 'linear'),
            ("tic8_Tmag", "dist_pc", '$T$ [mag]', '$d$ [pc]', 'linear', 'linear'),
            #('dec', 'ra', r'$\delta$ [deg]', r'$\alpha$ [deg]', 'linear', 'linear'),
            ('b', 'l', '$b$ [deg]', '$l$ [deg]', 'linear', 'linear'),
            ('tlc_mean_period', 'bp_rp', '$P$ [days]', '$G_{\mathrm{BP}}-G_{\mathrm{RP}}$ [mag]', 'log', 'linear'),
            #('rstar_sedfit', 'teff_sedfit', '$R_{\! \star}$ [$R_\odot$]', '$T_\mathrm{eff}$ [K]', 'linear', 'linear'),
            #('ruwe', 'bp_rp', 'RUWE', '$G_{\mathrm{BP}}-G_{\mathrm{RP}}$ [mag]', 'log', 'linear'),
            #('tlc_mean_period', 'Rcr_over_Rstar', '$P$ [hours]', '$R_{\mathrm{cr}}/R_{\! \star}$', 'log', 'linear'),
            ('Rcr_over_Rstar', 'tlc_mean_period',  '$R_{\mathrm{c}}/R_{\! \star}$', '$P$ [hours]', 'linear', 'log'),
            ('banyan_adopted_age', 'mass_parsec', 'Age [Myr]', '$M_{\! \star}$ [$M_\odot$]', 'log', 'linear'),
            #('rstar_sedfit', 'banyan_singleagefloat', '$R_{\! \star}$ [$R_\odot$]', 'Age [Myr]', 'linear', 'log'),
        ]
    elif plotsubset == 'sanitychecks':
        tuples = [
            ("M_G", "bp_rp", '$\mathrm{M}_{G}$ [mag]', '$G_{\mathrm{BP}}-G_{\mathrm{RP}}$ [mag]', 'linear', 'linear'),
            #("tic8_Tmag", "dist_pc", '$T$ [mag]', '$d$ [pc]', 'linear', 'linear'),
            #('dec', 'ra', r'$\delta$ [deg]', r'$\alpha$ [deg]', 'linear', 'linear'),
            #('b', 'l', '$b$ [deg]', '$l$ [deg]', 'linear', 'linear'),
            #('tlc_mean_period', 'bp_rp', '$P$ [days]', '$G_{\mathrm{BP}}-G_{\mathrm{RP}}$ [mag]', 'log', 'linear'),
            #('rstar_sedfit', 'teff_sedfit', '$R_{\! \star}$ [$R_\odot$]', '$T_\mathrm{eff}$ [K]', 'linear', 'linear'),
            #('tlc_mean_period', 'Rcr_over_Rstar', '$P$ [hours]', '$R_{\mathrm{cr}}/R_{\! \star}$', 'log', 'linear'),
            #('Rcr_over_Rstar', 'tlc_mean_period',  '$R_{\mathrm{cr}}/R_{\! \star}$', '$P$ [hours]', 'linear', 'log'),
            #('rstar_sedfit', 'banyan_singleagefloat', '$R_{\! \star}$ [$R_\odot$]', 'Age [Myr]', 'linear', 'log'),
            ('bp_rp', 'teff_sedfit', 'bp_rp', 'teff_sedfit', 'linear', 'linear'),
            ('banyan_adopted_age', 'mass_parsec', 'Age [Myr]', '$M_{\! \star}$ [$M_\odot$]', 'log', 'linear'),
            ('teff_sedfit', 'teff_parsec', 'teff_sedfit', 'teff_parsec', 'linear', 'linear'),
            ('rstar_sedfit', 'rstar_parsec', 'rstar_sedfit', 'rstar_parsec', 'linear', 'linear'),
            ('dist_metric_parsec', 'rstar_sedfit', 'dist_metric_parsec', 'rstar_sedfit', 'linear', 'linear'),
            ('dist_metric_parsec', 'mass_parsec', 'dist_metric_parsec', 'mass_parsec', 'linear', 'linear'),
            ('banyan_adopted_age', 'age_parsec', 'banyan_adopted_age', 'age_parsec', 'linear', 'linear'),
            ('dr3_ruwe', 'bp_rp', 'DR3 RUWE', '$G_{\mathrm{BP}}-G_{\mathrm{RP}}$ [mag]', 'log', 'linear'),
        ]
    else:
        raise NotImplementedError


    # make plot
    plt.close('all')
    set_style('clean')

    f = 1.6
    assert len(tuples) in [6,9]
    DO_NINECOLS = 0 if len(tuples) == 6 else 1

    if DO_NINECOLS:
        fig = plt.figure(figsize=(f*4,f*3))
        axd = fig.subplot_mosaic(
            """
            ABC
            DEF
            GHI
            """
        )
        axs = [axd[k] for k in 'A,B,C,D,E,F,G,H,I'.split(',')]
    else:
        fig = plt.figure(figsize=(f*3.5,f*2.5))
        axd = fig.subplot_mosaic(
            """
            ABC
            DEF
            """
        )
        axs = [axd[k] for k in 'A,B,C,D,E,F'.split(',')]

    for ix, ax in enumerate(axs):

        ykey, xkey, ylabel, xlabel, yscale, xscale = tuples[ix]

        # underplot the target sample distribution
        show_underplot = 0
        if ykey == 'M_G' and xkey == 'bp_rp':
            axkey = 'A'
            show_underplot = 1
            get_xval = lambda df: df[xkey]
            get_yval = lambda df: df[ykey]
        elif xkey == 'dist_pc' and ykey == 'tic8_Tmag':
            axkey = 'B'
            show_underplot = 1
            get_xval = lambda df: 1 / (df['parallax']*1e-3)
            get_yval = lambda df: df['TESSMAG']
        elif (xkey=='ra' and ykey=='dec') or (xkey=='l' and ykey=='b'):
            axkey = 'C'
            show_underplot = 1
            get_xval = lambda df: df[xkey]
            get_yval = lambda df: df[ykey]
        if show_underplot:
            ss = axd[axkey].get_subplotspec()
            axd[axkey].remove()
            import mpl_scatter_density # adds projection='scatter_density'
            axd[axkey] = fig.add_subplot(ss, projection='scatter_density')
            ax = axd[axkey]
            _d = underplot_cqvtargets(ax, get_xval, get_yval)
            if axkey == 'A':
                ax.text(0.97,0.97, 'CPVs', transform=ax.transAxes,
                        ha='right',va='top', color='C0', fontsize='small')
                ax.text(0.97,0.92, 'Candidates', transform=ax.transAxes,
                        ha='right',va='top', color='C0', fontsize='small',
                        alpha=0.5)
                ax.text(0.97,0.87, 'Searched TESS targets', transform=ax.transAxes,
                        ha='right',va='top', color='darkgray', fontsize='small',
                        alpha=0.9)
                ax.text(0.97,0.82, 'Unres. binaries', transform=ax.transAxes,
                        ha='right',va='top', color='orangered', fontsize='small',
                        alpha=1)


        if ykey == 'tlc_mean_period' and xkey == 'bp_rp':
            from gyrointerp.getters import get_Pleiades
            df_plei = get_Pleiades(overwrite=0)
            df_plei = df_plei[df_plei.flag_benchmark_period]
            ax.scatter(
                df_plei['dr2_bp_rp'], df_plei['Prot'], c='darkgray', s=3,
                linewidths=0, zorder=-1, alpha=0.7
            )
            ax.text(0.97,0.97, 'CPVs', transform=ax.transAxes,
                    ha='right',va='top', color='C0', fontsize='small')
            ax.text(0.97,0.92, 'Candidates', transform=ax.transAxes,
                    ha='right',va='top', color='lightblue', fontsize='small',
                    alpha=1)
            ax.text(0.97,0.87, 'Pleiades (R+16)', transform=ax.transAxes,
                    ha='right',va='top', color='darkgray', fontsize='small')
            ax.text(0.97,0.82, 'Unres. binaries', transform=ax.transAxes,
                    ha='right',va='top', color='orangered', fontsize='small',
                    alpha=1)


        fx = 1
        if xkey == 'tlc_mean_period' and ykey == 'Rcr_over_Rstar':
            fx = 24

        elif ykey == 'banyan_adopted_age':
            f0, f1, f2 = 1,1,1
            #np.random.seed(42)
            #ferr = np.random.normal(1, 0.02, size=len(df))
            #f0 = ferr*1.
            #print(f0)
            #np.random.seed(42)
            #ferr = np.random.normal(1, 0.02, size=len(maybe_df))
            #f1 = ferr*1.
            #print(f1)
            #np.random.seed(42)
            #ferr = np.random.normal(1, 0.02, size=len(binary_df))
            #f2 = ferr*1.
            #print(f2)

        else:
            f0 = 1
            f1 = 1
            f2 = 1

        ax.scatter(
            fx*df[xkey], f0*df[ykey], c='C0', s=6, linewidths=0, zorder=10
        )
        if showmaybe:
            ax.scatter(
                fx*maybe_df[xkey], f1*maybe_df[ykey], c='lightblue', s=6, linewidths=0, alpha=1, zorder=9
            )
        if emphruwe:
            ax.scatter(
                fx*binary_df[xkey], f2*binary_df[ykey], c='orangered', s=12, linewidths=0, alpha=1, zorder=8
            )

        ax.set_xscale(xscale)
        ax.set_yscale(yscale)

        if xkey == 'bp_rp':
            ax.set_xlim([1.4,4.5])

        if ykey == 'M_G':
            ylim = ax.get_ylim()[::-1]
            ax.set_ylim([15,5])
            assert max(df[ykey]) < 15
            assert min(df[ykey]) > 5
            ax.set_yticks([14,10,6])
            ax.set_yticklabels([14,10,6])

        if xkey == 'tlc_mean_period' and ykey == 'Rcr_over_Rstar':
            ax.set_xlabel("$P$ [hours]")
            ax.set_xlim([2, 48])
            assert max(df[xkey]) < 48
            assert max(maybe_df[xkey]) < 48

            # epic solution to the problem of "how do I get my mixed
            # major and minor labels working on a log plot?"
            minor_locator = FixedLocator([2,3,4,5,6,7,8,9,10,20,30,40])
            ax.xaxis.set_minor_locator(minor_locator)
            major_locator = FixedLocator([10])
            ax.xaxis.set_major_locator(major_locator)
            ax.set_xticklabels([10])
            def fn(x, pos):
                labels = [None,3,None,None,None,None,None,None,None,30,None]
                return labels[pos]
            formatter = FuncFormatter(fn)
            ax.xaxis.set_minor_formatter(formatter)

            ax.set_yticks([2,4,6])
            ax.set_yticklabels([2,4,6])
            minor_locator = FixedLocator([1,3,5,7])
            ax.yaxis.set_minor_locator(minor_locator)

        if ykey == 'tlc_mean_period' and xkey == 'bp_rp':
            ax.set_yticks([0.1, 1, 10])
            ax.set_yticklabels([0.1, 1, 10])

        if ykey == 'tic8_Tmag':
            ax.set_ylim([9,16])

        if xkey == 'dist_pc':
            ax.set_xlim([5, 155])
            assert max(df['dist_pc']) < 155
            assert max(maybe_df['dist_pc']) < 155

        if ykey == 'ruwe':
            ax.set_yticks([1, 10])
            ax.set_yticklabels([1, 10])

        if xkey == 'a_over_Rstar':
            minor_locator = FixedLocator([1,3,5,7])
            ax.xaxis.set_minor_locator(minor_locator)
            major_locator = FixedLocator([2,4,6])
            ax.xaxis.set_major_locator(major_locator)
            ax.set_xticklabels([2,4,6])

        if xkey == 'banyan_adopted_age':
            ax.set_xlim([0.9,250])
            ax.set_xticks([1,10,100])
            ax.set_xticklabels([1,10,100])

        if ykey == 'banyan_adopted_age':
            ax.set_ylim([0.9,250])
            ax.set_yticks([1,10,100])
            ax.set_yticklabels([1,10,100])

        if xkey == 'bp_rp':
            ax.set_xticks([2,3,4])
            ax.set_xticklabels([2,3,4])

        if ykey == 'rstar_sedfit':
            ax.set_yticks([0.2,0.6,1.0,1.4])
            ax.set_yticklabels([0.2,0.6,1.0,1.4])

        ax.set_xlabel(xlabel, labelpad=0.2)
        ax.set_ylabel(ylabel, labelpad=0.3)

        if (
            ('rstar' in xkey and 'rstar' in ykey)
            or
            ('age' in xkey and 'age' in ykey)
            or
            ('teff' in xkey and 'teff' in ykey)
        ):
            xlo, xhi = ax.get_xlim()
            ylo, yhi = ax.get_ylim()
            ax.plot([0, max([xhi, yhi])], [0, max([xhi, yhi])], c='k', lw=0.5, zorder=-100, alpha=0.5)
            ax.set_xlim([xlo, xhi])
            ax.set_ylim([ylo, yhi])


    # set naming options
    s = ''
    if showmaybe:
        s += '_showmaybe'
    if isinstance(plotsubset, str):
        s += f'_{plotsubset}'
    if emphruwe:
        s += '_emphruwe'

    bn = 'catalogscatter'
    outpath = join(outdir, f'{bn}{s}.png')

    fig.tight_layout(w_pad=0.2, h_pad=0.2)

    savefig(fig, outpath, dpi=400)


def underplot_cqvtargets(ax, get_xval, get_yval):

    from complexrotators.getters import get_cqv_search_sample
    df_bkgd = get_cqv_search_sample()

    from matplotlib.colors import LinearSegmentedColormap
    # "Viridis-like" colormap with white background
    white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
        (0, '#ffffff'),
        (1e-20, '#440053'),
        (0.2, '#404388'),
        (0.4, '#2a788e'),
        (0.6, '#21a784'),
        (0.8, '#78d151'),
        (1, '#fde624'),
    ], N=256)

    _x = get_xval(df_bkgd)
    _y = get_yval(df_bkgd)
    s = np.isfinite(_x) & np.isfinite(_y)

    # add log stretch...
    from astropy.visualization import LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    vmin, vmax = 0.1, 10000
    vmin, vmax = 0.1, 5000

    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())

    cmap = "Greys"
    density = ax.scatter_density(_x[s], _y[s], cmap=cmap, norm=norm)

    return density


def plot_magnetic_bstar_comparison(outdir, showtitles=1, titlefontsize=3.75,
                                   selfn='hd37776'):

    # mostmassiveCQV,  HD 37776/landstreet/V901 Ori
    if selfn == 'hd37776':
        ticids = ['11400909', '405754448']
        optionalid = ['HD 37776', None]
        showsectors = [6, 38]
        texts = [7, 0.8] # masses
    # HD 64740, in VelaOB2
    elif selfn == 'hd64740':
        ticids = ['268971806', '201789285']
        optionalid = ['HD 64740', None]
        showsectors = [34, 3]
        texts = [7, 0.1] # masses
    elif selfn == 'both':
        ticids = ['405754448', '11400909', '201789285', '268971806']
        optionalid = [None, 'HD 37776', None, 'HD 64740']
        showsectors = [38, 6, 3, 34]
        texts = [0.8, 7, 0.1, 7] # masses
    elif selfn == 'all':
        ticids = [
            '11286209', '424048289', # simple B star field topplogies
            '11400909', '268971806', # complex B star field topologies
            '405754448', '201789285', # M dwarf matches
        ]
        optionalid = [
            'σ Ori E', 'HD 345439',
            #None, 'HD 37776',
            #None, 'HD 64740'
            'HD 37776', 'HD 64740',
            None, None
        ]
        showsectors = [
            32, 41,
            #38, 6,
            #3, 34
            6, 34,
            38, 3
        ]
        # sigmaOriE from Townsend2013
        texts = [
            8, 7,
            7, 7,
            0.8, 0.1
        ] # masses
    elif selfn == 'simpleB':
        ticids = [
            '11286209', '424048289' # simple B star field topplogies
        ]
        optionalid = [
            'σ Ori E', 'HD 345439',
        ]
        showsectors = [
            32, 41,
        ]
        # sigmaOriE from Townsend2013
        texts = [
            8, 7,
        ] # masses
    elif selfn == 'complexB':
        ticids = [
            '11400909', '268971806', # complex B star field topologies
        ]
        optionalid = [
            'HD 37776', 'HD 64740',
        ]
        showsectors = [
            6, 34,
        ]
        # sigmaOriE from Townsend2013
        texts = [
            7, 7,
        ] # masses
    elif selfn == 'complexM':
        ticids = [
            '405754448', '201789285', # M dwarf matches
        ]
        optionalid = [
            None, None
        ]
        showsectors = [
            38, 3
        ]
        # sigmaOriE from Townsend2013
        texts = [
            0.8, 0.1
        ] # masses


    # sigma-ori
    elif selfn == 'sigmaoriE':
        msg = (
            "there are not any good matches..."
            "because i would have labelled sigma Ori E an RS CVn"
        )
        raise NotImplementedError(msg)
    else:
        raise NotImplementedError

    # prepare to get lc data
    from complexrotators.lcprocessing import (
        cpv_periodsearch, prepare_cpv_light_curve
    )

    #
    # instantiate plot
    #
    set_style("science")
    if len(ticids) == 2:
        factor = 1.4
        fig, axs = plt.subplots(figsize=(factor*2,factor*1.37), ncols=2)
    elif len(ticids) == 4:
        fig, axs = plt.subplots(figsize=(4,1.5), ncols=4)
    elif len(ticids) == 6:
        fig, axs = plt.subplots(figsize=(6,1.5), ncols=6)

    ix = 0

    for ticid, sector, ax, opt in zip(ticids, showsectors, axs, optionalid):

        lcdir = '/Users/luke/.lightkurve/cache/mastDownload/TESS'
        lcpaths = glob(join(
            lcdir, f"tess*s{str(sector).zfill(4)}*{ticid}*",
            f"*{ticid}*_lc.fits")
        )
        if not len(lcpaths) == 1:
            print('bad lcpaths')
            print(ticid, sector)
            print(lcpaths)
            import IPython; IPython.embed()
            assert 0

        lcpath = lcpaths[0]

        cachedir = outdir

        # get the relevant light curve data
        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid) = prepare_cpv_light_curve(lcpath, cachedir)

        if ticid == '424048289':
            mask = (
                (x_obs < 2423) |
                ( (x_obs > 2433) & (x_obs < 2437) ) |
                (x_obs > 2444) |
                ( (x_obs > 2430) & (x_obs < 2435) )
            )
            sel = ~mask
            x_obs = x_obs[sel]
            y_flat = y_flat[sel]
            #plt.close('all')
            #plt.figure(figsize=(12,4))
            #plt.scatter(x_obs, y_flat, s=1)
            #plt.savefig('temp.png')
            #assert 0

        if ticid == '11400909':
            sel = x_obs < 1487
            x_obs = x_obs[sel]
            y_flat = y_flat[sel]

        # get period, t0, and periodogram (PDM or LombScargle)
        d = cpv_periodsearch(
            x_obs, y_flat, starid, cachedir, t0='binmin', periodogram_method='pdm'
        )
        if ticid == '11400909':
            d['period'] *= 0.5

        bd = time_bin_magseries(d['times'], d['fluxs'], binsize=1200, minbinelems=1)

        ylim = get_ylimguess(1e2*(bd['binnedmags']-np.nanmean(bd['binnedmags'])))

        if showtitles:
            if isinstance(opt, str):
                titlestr = f'{opt}, S{sector}, {d["period"]*24:.1f}h'
            else:
                titlestr = f'TIC {ticid}, S{sector}, {d["period"]*24:.1f}h'
        else:
            titlestr = None

        binsize_phase = 1/300

        BINMS=1.0
        alpha0=0.15

        yoffset = 0
        if ticid == '405754448':
            yoffset = 0.25

        titlefontsize = 'xx-small'
        plot_phased_light_curve(
            d['times'], d['fluxs'], d['t0'], d['period'], None, ylim=ylim,
            xlim=[-0.6,0.6], binsize_phase=binsize_phase, BINMS=BINMS, titlestr=titlestr,
            showtext=None, showtitle=False, figsize=None, c0='darkgray',
            alpha0=alpha0, c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
            fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
            showxticklabels=False, titlefontsize=titlefontsize,
            titlepad=0.05, yoffset=yoffset, rasterized=rasterized
        )

        if not showtitles:
            txt = f'{d["period"]*24:.1f}h'
            tform = ax.transAxes
            props = dict(boxstyle='square', facecolor='white', alpha=0.7, pad=0.15,
                         linewidth=0)
            fontsize = 'xx-small'
            ax.text(
                0.97, 0.05, txt, transform=tform, ha='right',
                va='bottom', color='k', fontsize=fontsize, bbox=props
            )

        if isinstance(texts, list):
            text = texts[ix]
            #txt = r'$\approx \!$' + f'{text}' + '$\,M_\odot$'
            txt = f'{text}' + '$\,M_\odot$'
            tform = ax.transAxes
            props = dict(boxstyle='square', facecolor='white', alpha=1, pad=0.15,
                         linewidth=0)
            fontsize = 'xx-small'
            ax.text(
                0.98, 0.05, txt, transform=tform, ha='right', va='bottom',
                color='k', fontsize=fontsize, bbox=props
            )

        ax.set_xticks([-0.5,0,0.5])

        ylow, yhigh = int(np.ceil(ylim[0]))+1, int(np.floor(ylim[1]))-1
        if np.diff([np.abs(ylow), yhigh]) <= 2:
            ylowabs = np.abs(ylow)
            yhighabs = np.abs(yhigh)
            ylow = -np.min([ylowabs, yhighabs])
            yhigh = np.min([ylowabs, yhighabs])
            if ylow == yhigh == 0:
                # 368129164 led to this
                ylow = -1
                yhigh = 1
            if yhigh >= 10:
                ylow = -9
                yhigh = 9

        ax.set_yticks([ylow, 0, yhigh])
        ax.set_yticklabels([ylow, 0, yhigh])

        id_yticks = {
            "405754448": [-1,0,1],
        }
        if str(ticid) in list(id_yticks.keys()):
            yticks = id_yticks[str(ticid)]
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks)
            if yticks[2] >= 9:
                diff = 4
            elif yticks[2] > 5:
                diff = 2.5
            elif yticks[2] > 2:
                diff = 1.5
            else:
                diff = 1.
            ax.set_ylim([yticks[0]-diff, yticks[2]+diff])

            if str(ticid) == '397791443':
                ax.set_ylim([-17,13])

        if ticid == '405754448':
            ax.set_ylim([-1.5, 1.5])
            ax.set_yticks([-1, 0, 1])
            ax.set_yticklabels([-1, 0, 1])
        if ticid == '11400909': # hd37776
            ax.set_ylim([-1.5, 1.5])
            ax.set_yticks([-1, 0, 1])
            ax.set_yticklabels([-1, 0, 1])
        if ticid == '201789285':
            ax.set_ylim([-7, 7])
            ax.set_yticks([-5, 0, 5])
            ax.set_yticklabels([-5, 0, 5])
        if ticid == '268971806': # hd64740
            ax.set_ylim([-0.2, 0.2])
            ax.set_yticks([-0.1, 0, 0.1])
            ax.set_yticklabels([-0.1, 0, 0.1])
        if ticid == '11286209': # sigma ori E
            ax.set_ylim([-6, 4])
            ax.set_yticks([-5, 0, 3])
            ax.set_yticklabels([-5, 0, 3])
        if ticid == '424048289':
            ax.set_ylim([-6, 6])
            ax.set_yticks([-5, 0, 5])
            ax.set_yticklabels([-5, 0, 5])

        labelsize = 'xx-small'

        ax.minorticks_off()

        ax.tick_params(axis='both', which='major', labelsize=labelsize,
                       pad=1.5)

        ix += 1

    for ax in axs:
        ax.set_xticklabels(['-0.5','0','0.5'])

    fs = 'small'
    fig.text(0.5, 0, r"Phase, φ", fontsize=fs, va='bottom', ha='center')
    fig.text(-0.01,0.5, "$\Delta$ Flux [%]", va='center', rotation=90, fontsize=fs)

    # set naming options
    s = ''
    if showtitles:
        s += '_showtitles'
    if isinstance(selfn, str):
        s += f'_{selfn}'

    # height/width
    fig.tight_layout(h_pad=0.2, w_pad=0.2)

    outpath = join(outdir, f'magnetic_bstar_comparison{s}.png')
    savefig(fig, outpath, dpi=400)


def plot_lineevolnpanel(outdir, starid=None, jstr=None):
    """
    starid: e.g. 'TIC402980664'
    jstr: e.g. 'j531',
    """

    # get data
    lines = [
        #'Hη',
        'Hζ',
        'Ca[K]',
        'Ca[H] + Hε',
        'CaI',
        'Hγ',
        'He',
        'Na D1',
        'Hα'
    ]
    globs = [
        #'*bj*order04*',
        '*bj*order05*',
        '*bj*order07*',
        '*bj*order07*',
        '*bj*order13*',
        '*bj*order15*',
        '*rj*order10*',
        '*rj*order11*',
        '*ij*order00*'
    ]
    if starid == 'TIC402980664':
        deltawav = 5
    elif starid == 'TIC141146667':
        deltawav = 7
    xlims = [
        #[3835.397-1.5*deltawav, 3835.397+1.5*deltawav], # Hη (9->2)
        [3889.096-1.5*deltawav, 3889.096+1.5*deltawav], # Hζ (8->2)
        [3933.66-deltawav, 3933.66+deltawav], # ca k
        [3968.47-deltawav, 3968.47+deltawav], # ca h
        [4226-deltawav, 4226+deltawav], # CaI
        [4340.05-deltawav, 4340.05+deltawav], # Hγ
        [5875.618-deltawav, 5875.618+deltawav], #He
        [5895.92-deltawav, 5895.92+deltawav], # Na D1
        [6562.8-1.5*deltawav, 6562.8+1.5*deltawav], # halpha
    ]

    #
    # get all wavelengths & fluxes over all chips...
    #

    from cdips_followup.spectools import read_hires
    specdir = f'/Users/luke/Dropbox/proj/cpv/data/spectra/HIRES/{starid}_RDX'
    bpaths = np.sort(glob(join(specdir, f'b{jstr}*.fits')))
    ipaths = np.sort(glob(join(specdir, f'i{jstr}*.fits')))
    rpaths = np.sort(glob(join(specdir, f'r{jstr}*.fits')))

    # make plot
    plt.close('all')
    set_style('clean')

    toprowkeys = [x for x in 'abcdefgh']
    bottomrowkeys = [x for x in 'ijklmnop']

    #fig, ax = plt.subplots(figsize=(4,3))
    fig = plt.figure(figsize=(11*1.3,3*1.3))
    axd = fig.subplot_mosaic(
        """
        abcdefgh
        ijklmnop
        """,
        gridspec_kw={
            #"width_ratios": [1, 1, 1, 1]
            "height_ratios": [1, 3]
        },
    )

    from scipy.ndimage import gaussian_filter1d
    from specutils.spectra import Spectrum1D, SpectralRegion
    from specutils.fitting import fit_generic_continuum
    import warnings
    from rudolf.plotting import multiline

    # iterate over each line / spectral window
    for ix, (l,g,xlim) in enumerate(zip(lines, globs, xlims)):

        toprowkey = toprowkeys[ix]
        bottomrowkey = bottomrowkeys[ix]

        if 'bj' in g:
            specpaths = bpaths
        elif 'ij' in g:
            specpaths = ipaths
        elif 'rj' in g:
            specpaths = rpaths

        # iterate over times...
        mjds, flxs, wavs = [],[],[]
        for specpath in specpaths:

            # header time info
            hdul = fits.open(specpath)
            dateobs = hdul[0].header['DATE-OBS']
            exptime = hdul[0].header['EXPTIME']
            mjd = hdul[0].header['MJD']
            utc = hdul[0].header['UTC'][:5]
            timestr = f"{dateobs} {utc}UT"
            hdul.close()

            # flux info
            flx_2d, wav_2d = read_hires(
                specpath, is_registered=0, return_err=0
            )

            # there are a few cases, like CaHK, for which multiple orders exist
            order = int(g.split('order')[-1].rstrip("*"))
            flx, wav = flx_2d[order, :], wav_2d[order, :]

            #flx, wav = get_flx_wav_given_2d_and_target(
            #    flx_2d, wav_2d, xlim[1]-deltawav
            #)

            # "normalize"
            #flx /= np.nanmedian(flx)
            if l != 'He':
                flx /= np.nanmedian(flx)
            else:
                _sel = (wav < 5881) & (wav > 5879)
                normfactor = np.nanmedian(flx[_sel])
                flx /= normfactor

            mjds.append(mjd)
            flxs.append(flx)
            wavs.append(wav)

            ## range over which to continuum normalize...
            #sel = ((xlim[0]-20) < wav) & (wav < xlim[1]+20)

            #swav = wav[sel]
            #sflx = flx[sel]

            # # continuum normalize
            # spec = Spectrum1D(spectral_axis=swav*u.AA,
            #                   flux=sflx*u.dimensionless_unscaled)
            # with warnings.catch_warnings():  # Ignore warnings
            #     warnings.simplefilter('ignore')
            #     g1_fit = fit_generic_continuum(spectrum)
            # sflx_continuum_fitted = g1_fit(sflx*u.AA)
            # spec_normalized = spec / sflx_continuum_fitted
            # norm_sflx = spec_normalized.flux
            # spectrum = Spectrum1D(flux=y*u.dimensionless, spectral_axis=x*u.um)

        mjds = np.array(mjds).astype(float)
        flxs = np.array(flxs)
        wavs = np.array(wavs)

        # does the wavelength solution vary over time?  it should not.  this
        # assertion statement verifies that.
        assert np.diff(wavs, axis=0).sum() == 0

        # average over time to get the median line profile
        # compute direct subtraction residual
        fn = lambda x: gaussian_filter1d(x, sigma=4)

        flx_median = np.nanmedian(flxs, axis=0)
        smooth_flx_median = fn(flx_median)

        smooth_flxs = np.array([
            fn(flxs[ix, :]) for ix in range(len(flxs))
        ])

        smooth_diff_flxs = smooth_flxs - smooth_flx_median[None, :]

        # plot average line profile over all spectra
        sel = ((xlim[0]-1.3*deltawav) < wav) & (wav < xlim[1]+1.3*deltawav)
        swav = wav[sel]
        sflx = smooth_flx_median[sel]

        axd[toprowkey].plot(
            swav, sflx, c='k', zorder=3, lw=0.2
        )

        # plot individual line profiles at each time (minus the average line
        # profile)
        lc = multiline(
            wavs, smooth_diff_flxs, 24*(nparr(mjds)-np.min(mjds)),
            #cmap='Spectral',
            #cmap='winter',
            cmap='viridis',
            ax=axd[bottomrowkey], lw=0.5
        )

        axd[toprowkey].set_title(l)
        axd[toprowkey].set_xticklabels([])

        axd[toprowkey].set_xlim(xlim)
        axd[bottomrowkey].set_xlim(xlim)

        axd[bottomrowkey].xaxis.set_major_locator(MaxNLocator(5))

        if starid == 'TIC402980664':
            if l == 'Hη':
                axd[bottomrowkeys[ix]].set_ylim([-2.5, 8])
            elif l == 'CaI':
                axd[bottomrowkeys[ix]].set_ylim([-0.5, 2])
            elif l == 'He':
                axd[bottomrowkeys[ix]].set_ylim([-0.3, 1.4])
            elif l == 'Na D1':
                axd[bottomrowkeys[ix]].set_ylim([-0.3, 0.9])
            elif l == 'Hζ':
                axd[bottomrowkeys[ix]].set_ylim([-3, 11])
            #if ix == 2:
            #    axd[str(ix+3)].set_yticks([-0.1,0,0.1])
        if starid == 'TIC141146667':
            axd[bottomrowkeys[ix]].set_ylim([-10, 10])
            if l == 'Hα':
                axd[bottomrowkeys[ix]].set_ylim([-2, 2])
            elif l == 'Hη':
                axd[bottomrowkeys[ix]].set_ylim([-2.5, 8])
            elif l == 'CaI':
                axd[bottomrowkeys[ix]].set_ylim([-2, 4])
            elif l == 'He':
                axd[bottomrowkeys[ix]].set_ylim([-1, 1])
            elif l == 'Na D1':
                axd[bottomrowkeys[ix]].set_ylim([-1, 1])
            elif l == 'Hζ':
                axd[bottomrowkeys[ix]].set_ylim([-10, 10])

        expids = [os.path.basename(s).rstrip('.fits').lstrip('b')
                  for s in specpaths]
        outdict = {
            'wav': wav,
            'flxs': flxs,
            'flx_median': flx_median,
            'smooth_flx_median': smooth_flx_median,
            'smooth_diff_flxs': smooth_diff_flxs,
            'mjds': mjds,
            'specpaths': specpaths,
            'expids': expids
        }
        cachepath = join(
            outdir, f'spec_cache_{g.replace("*","_").lstrip("_").rstrip("_")}.pkl'
        )
        with open(cachepath, "wb") as f:
            pickle.dump(outdict, f)
        print(f'Wrote {cachepath}')


    axd[toprowkeys[0]].set_ylabel('Median line profile')
    axd[bottomrowkeys[0]].set_ylabel('Observed - Median line profile')

    axins1 = inset_axes(axd[bottomrowkeys[-1]], width="25%", height="5%",
                        loc='lower right', borderpad=1.5)
    cb = fig.colorbar(lc, cax=axins1, orientation="horizontal")
    cb.ax.tick_params(labelsize='xx-small')
    cb.ax.set_title('Time [hours]', fontsize='xx-small')
    # axd['b'].text(0.725,0.955, '$t$ [days]',
    #         transform=axd['b'].transaxes,
    #         ha='right',va='top', color='k', fontsize='xx-small')
    cb.ax.tick_params(size=0, which='both') # remove the ticks
    axins1.xaxis.set_ticks_position("bottom")

    #fig.text(-0.01,0.5, 'Flux [order-normalized median]', va='center',
    #         rotation=90)
    fig.text(0.5,-0.01, r'Wavelength [$\AA$]', va='center', ha='center',
             rotation=0)

    fig.tight_layout()

    # set naming options
    s = ''

    outpath = os.path.join(outdir, f'{starid}_{jstr}_lineevolnpanel.png')
    savefig(fig, outpath, dpi=400)


def plot_movie_phase_timegroups(
    outdir,
    ticid=None,
    lc_cadences='2min',
    manual_period=None,
    t0='binmin',
    ylim=None,
    binsize_phase=0.005,
    xlim=[-0.6,0.6],
    yoffset=5,
    showtitle=1,
    figsize_y=7,
    model_id=None,
    rasterized=False,
    N_cyclestobin=3,
    sector_range=None,
    style='science',
    arial_font=0
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
    _times, _fluxs, _t0s, _periods, _titlestrs, _sectorstrs = [],[],[],[],[],[]

    for lc in lclist:

        if sector_range is not None:
            sector = lc.sector
            ok_sector_list = list(sector_range)
            if sector not in ok_sector_list:
                continue

        (time, flux, qual, x_obs, y_obs, y_flat,
         y_trend, x_trend, cadence_sec, sector,
         starid) = prepare_given_lightkurve_lc(lc, ticid, outdir)

        # get t0, period, lsp
        if not isinstance(t0, float) and isinstance(manual_period, float):
            d = cpv_periodsearch(x_obs, y_flat, starid, outdir, t0=t0)
        else:
            d = {'times': x_obs, 'fluxs': y_flat,
                 't0': t0, 'period': manual_period}

        _t0 = d['t0']
        if isinstance(t0, float):
            _t0 = t0
        period = d['period']
        if isinstance(manual_period, float):
            period = manual_period
        titlestr = starid.replace('_',' ')

        _times.append(d['times'])
        _fluxs.append(d['fluxs'])
        _t0s.append(_t0)
        _periods.append(period)
        _titlestrs.append(titlestr)
        _sectorstrs.append(np.repeat(sector, len(d['times'])))

    # merge lightcurve data, and split before making the plot.
    times = np.hstack(_times)
    fluxs = np.hstack(_fluxs)
    t0s = np.hstack(_t0s)
    periods = np.hstack(_periods)
    titlestrs = np.hstack(_titlestrs)
    sectorstrs = np.hstack(_sectorstrs)

    if isinstance(model_id, str):
        # for now, just lp12-502
        assert ticid == "TIC_402980664"
        manual_csvpath = f'/Users/luke/Dropbox/proj/cpv/results/4029_mask/lc_lsresid_{model_id}.csv'
        df = pd.read_csv(manual_csvpath)
        times = np.array(df.time)
        # residual flux from subtracting the model given in model_id
        fluxs = np.array(df.r_flux)

    #
    # determine break times and groups
    #
    plot_period = np.nanmean(_periods)
    plot_t0 = t0s[0]

    min_cycle = 0 # NOTE: this might cut depending on t0
    max_cycle = int(np.ceil((np.nanmax(times) - plot_t0)/plot_period))
    break_times = [plot_t0 + ix*plot_period for ix in
                   range(min_cycle, max_cycle, N_cyclestobin)]

    # big time groups, separated by at least like a week.
    # (might use these if you want to generate gap frames??)
    from astrobase.lcmath import find_lc_timegroups
    ngroups, groups = find_lc_timegroups(times, mingap=7*24/24)

    # Make plots
    time_index = 0
    frames_per_second = 3
    N_gap_seconds = 1.
    N_frames_for_gaps = int(np.ceil(frames_per_second * N_gap_seconds))
    timegap_counter = 0
    last_e_end = 0

    for t_start, t_stop in zip(break_times[:-1], break_times[1:]):

        sel = (times >= t_start) & (times <= t_stop)
        N_times = len(times[sel])

        e_start = int(np.floor((t_start - plot_t0)/plot_period))
        e_end = int(np.floor((t_stop - plot_t0)/plot_period))
        if last_e_end == 0:
            last_e_end == e_end

        # require at least ~one cycle, unless we are at the beginning of a
        # time gap.  in those instances, plot "Data Gap" frames for
        # N_gap_seconds.
        VERBOSE = 0
        if VERBOSE:
            print(time_index, t_start, t_stop, N_times, last_e_end, e_end, timegap_counter)

        FLAG_TIMEGAP = 0
        if N_times <= 0.9*plot_period*24*60/2:
            if e_end - last_e_end > 30 and timegap_counter < N_frames_for_gaps:
                # in this instance, there is a time gap.
                FLAG_TIMEGAP = 1
            else:
                continue
        if timegap_counter == N_frames_for_gaps:
            timegap_counter = 0

        iso_t0 = Time(t_start+2457000, format='jd').isot[:10]
        iso_t1 = Time(t_stop+2457000, format='jd').isot[:10]

        gtime = times[sel]
        gflux = fluxs[sel]
        gsectors = sectorstrs[sel]
        assert len(np.unique(gsectors)) <= 2
        try:
            _sector = gsectors[0]
        except IndexError:
            pass

        plt.close('all')
        set_style(style)
        if arial_font:
            rcParams['font.family'] = 'Arial'
        factor=1.
        fig, ax = plt.subplots(figsize=(factor*3, factor*3))

        txt0 = f"{iso_t0}"+"$\,$-$\,$"+f"{iso_t1}"
        txt1 = f"Cycle {e_start}"+"$\,$-$\,$"+f"{e_end}, Sector {_sector}"
        txt = txt0 + '\n' + txt1

        c1 = 'k' if 'wob' not in style else 'white'
        if not FLAG_TIMEGAP:
            plot_phased_light_curve(
                gtime, gflux, plot_t0, plot_period, None,
                fig=fig, ax=ax,
                binsize_phase=binsize_phase,
                xlim=xlim,
                #showtext=txt,
                c0='darkgray',
                c1=c1,
                titlestr=txt,
                titlepad=0.1,
                showtext=False,
                savethefigure=False,
                titlefontsize='xx-small',
                rasterized=rasterized
            )
        if isinstance(ylim, (list, tuple)):
            ax.set_ylim(ylim)
        if isinstance(xlim, (list, tuple)):
            ax.set_xlim(xlim)

        if FLAG_TIMEGAP:
            ax.text(
                0.5, 0.5, 'Gap', ha='center', va='center',
                fontsize='x-large', transform=ax.transAxes
            )
            ax.set_title(txt, fontsize='xx-small', pad=0.1, color='k')
            ax.set_xticks([-0.5, -0.25, 0, 0.25, 0.5])
            ax.set_xticklabels([-0.5, -0.25, 0, 0.25, 0.5])

        ax.set_ylabel(r"$\Delta$ Flux [%]", fontsize='large')
        ax.set_xlabel(r"Phase, φ", fontsize='large')

        format_ax(ax)
        fig.tight_layout()

        s = ''
        if isinstance(model_id, str):
            s += f'_{model_id}'
        if rasterized:
            s += "_rasterized"
        if 'wob' in style:
            s += '_wob'

        tstr = str(time_index).zfill(4)

        outpath = join(
            outdir,
            f"{ticid}_{tstr}_{lc_cadences}_phase_timegroups{s}.png"
        )

        fig.savefig(outpath, bbox_inches='tight', dpi=450)
        print(f"saved {outpath}")

        time_index += 1
        if not FLAG_TIMEGAP:
            last_e_end = e_end

        if FLAG_TIMEGAP:
            if timegap_counter < N_frames_for_gaps:
                timegap_counter += 1
            else:
                last_e_end = e_end
                timegap_counter = 0


def plot_movie_specriver(
    outdir,
    ticid=None,
    linestr='Hα',
    lc_cadences='2min',
    manual_period=None,
    t0='binmin',
    ylim=None,
    binsize_phase=0.005,
    xlim=[-0.6,0.6],
    showtitle=1,
    figsize_y=7,
    rasterized=False,
    sector: int=None,
    style='science',
    arial_font=0,
    lamylim=None,
    cb_ticks=[1,2,3]
    ):
    """
    As in plot_phase
    """

    ###############
    # get LC data #
    ###############
    assert isinstance(sector, int)
    sector_range = [sector]

    lclist = _get_cpv_lclist(lc_cadences, "TIC "+ticid)

    if len(lclist) == 0:
        print(f'WRN! Did not find light curves for {ticid}. Escaping.')
        return 0

    # for each light curve (sector / cadence specific), detrend if needed, get
    # the best period.
    _times, _fluxs, _t0s, _periods, _titlestrs, _sectorstrs = [],[],[],[],[],[]

    for lc in lclist:

        if sector_range is not None:
            sector = lc.sector
            ok_sector_list = list(sector_range)
            if sector not in ok_sector_list:
                continue

        (time, flux, qual, x_obs, y_obs, y_flat,
         y_trend, x_trend, cadence_sec, sector,
         starid) = prepare_given_lightkurve_lc(lc, ticid, outdir)

        # get t0, period, lsp
        if not isinstance(t0, float) and isinstance(manual_period, float):
            d = cpv_periodsearch(x_obs, y_flat, starid, outdir, t0=t0)
        else:
            d = {'times': x_obs, 'fluxs': y_flat,
                 't0': t0, 'period': manual_period}

        _t0 = d['t0']
        if isinstance(t0, float):
            _t0 = t0
        period = d['period']
        if isinstance(manual_period, float):
            period = manual_period
        titlestr = starid.replace('_',' ')

        _times.append(d['times'])
        _fluxs.append(d['fluxs'])
        _t0s.append(_t0)
        _periods.append(period)
        _titlestrs.append(titlestr)
        _sectorstrs.append(np.repeat(sector, len(d['times'])))

    # merge lightcurve data, and split before making the plot.
    times = np.hstack(_times)
    fluxs = np.hstack(_fluxs)
    t0s = np.hstack(_t0s)
    periods = np.hstack(_periods)
    titlestrs = np.hstack(_titlestrs)
    sectorstrs = np.hstack(_sectorstrs)

    #################
    # get spec data #
    #################
    from complexrotators.getters import get_specriver_data
    specpaths, spectimes, xvals, yvals = get_specriver_data(
        ticid, linestr
    )

    t0 = t0s[0]
    period = periods[0]

    _pd = phase_magseries(spectimes, np.ones(len(spectimes)), period, t0,
                          wrap=0, sort=False)
    specphases = _pd['phase']
    specphases[specphases > 0.5] -= 1

    flux_arr = np.zeros(
        (len(xvals[0]), len(spectimes))
    )
    for ix, yval in enumerate(yvals):
        flux_arr[:, ix] = yval


    ##########################################
    # Begin the plot

    for time_index, (spectime, specphase, specpath, xval, yval) in enumerate(zip(
        spectimes, specphases, specpaths, xvals, yvals
    )):

        plt.close("all")
        set_style(style)
        if arial_font:
            rcParams['font.family'] = 'Arial'

        fig = plt.figure(figsize=(8,3))
        axd = fig.subplot_mosaic(
            """
            ABC
            """#,
            #gridspec_kw={
            #    "width_ratios": [1, 1, 1, 1]
            #}
        )

        ##########################################
        # flux vs phase
        ax = axd['A']
        txt = ''
        c0 = 'darkgray'
        c1 = 'k' if 'wob' not in style else 'white'
        plot_phased_light_curve(
            times, fluxs, t0, period, None, fig=fig, ax=ax, titlestr=None,
            binsize_phase=binsize_phase, xlim=xlim, yoffset=0, showtext=txt,
            savethefigure=False, dy=0, rasterized=rasterized, c0=c0, c1=c1
        )
        if isinstance(ylim, (list, tuple)):
            ax.set_ylim(ylim)
        if isinstance(xlim, (list, tuple)):
            ax.set_xlim(xlim)

        ylim = ax.get_ylim()
        ax.vlines(specphase, ylim[0], ylim[1], colors='darkgray', alpha=0.5,
                  linestyles='--', zorder=-10, linewidths=0.5)
        ax.set_ylim(ylim)

        ax.set_ylabel(r"$\Delta$ Flux [%]", fontsize='large')
        ax.set_xlabel(r"Phase, φ", fontsize='large')
        format_ax(ax)

        ##########################################
        # flux vs wavelength
        ax = axd['B']

        c = 'k' if 'wob' not in style else 'white'
        ax.plot(
            xval, yval, c=c, lw=0.5
        )
        txt = f't={24*(spectime-min(spectimes)):.1f}hr, φ={specphase:.2f}'
        ax.text(
            0.96, 0.96, txt, ha='right', va='top', transform=ax.transAxes
        )
        txt = linestr
        ax.text(
            0.04, 0.96, txt, ha='left', va='top', transform=ax.transAxes
        )

        assert isinstance(lamylim, (list, tuple))
        ax.set_ylim(lamylim)

        ax.set_ylabel("$f_\lambda$", fontsize='large')
        ax.set_xlabel(r"Δv [km/s]", fontsize='large')

        ##########################################
        # specriver: phase vs wavelength, color by flux
        ax = axd['C']

        cmap = 'YlGnBu'
        cmap = 'Greys_r'
        vmin = lamylim[0]
        vmax = lamylim[1]
        c = ax.pcolor(xval,
                      24*(spectimes-min(spectimes)),
                      flux_arr.T,
                      cmap=cmap,
                      norm=colors.LogNorm(
                           vmin=0.9, vmax=vmax
                      ),
                      shading='auto', rasterized=True)

        ax.set_ylabel("Time [hr]", fontsize='large')
        ax.set_xlabel(r"Δv [km/s]", fontsize='large')

        xmin, xmax = ax.get_xlim()
        ax.hlines(24*(spectime-min(spectimes)), xmin, xmax,
                  colors='darkgray', alpha=0.9, linestyles='--', zorder=2,
                  linewidths=0.5)
        ax.set_xlim((xmin, xmax))

        # sick inset colorbar
        x0,y0,dx,dy = 1.02, -0.09, 0.3, 0.02
        axins1 = inset_axes(ax, width="100%", height="100%",
                            bbox_to_anchor=(x0,y0,dx,dy),
                            loc='lower right',
                            bbox_transform=ax.transAxes)
        cb = fig.colorbar(c, cax=axins1, orientation="horizontal",
                          extend="both")

        cb.set_label("$f_\lambda$", rotation=0, labelpad=3)
        if cb_ticks is not None:
            cb.set_ticks(cb_ticks)
            cb.set_ticklabels(cb_ticks)
        else:
            cb.ax.yaxis.set_tick_params(left=False, labelleft=False)
            cb.ax.xaxis.set_tick_params(left=False, labelleft=False)
            cb.ax.xaxis.set_tick_params(bottom=False, labelbottom=False)

            # NOTE doesn't work...
            # cb.set_ticks([])
            # tick_labels = cb.ax.get_yticklabels()
            # #plt.setp(tick_labels, visible=False)
            # #cb.ax.tick_params(color='k')
            # ## Get the current tick locations and labels
            # #tick_locs = [1,vmax]
            # #tick_labels = [1, int(vmax)]
            # ## Set the new tick locations and labels
            # #cb.set_ticks(tick_locs)
            # cb.set_ticklabels(tick_labels, fontsize='xx-small')


        fig.tight_layout()

        s = ''
        if rasterized:
            s += "_rasterized"
        if 'wob' in style:
            s += '_wob'

        tstr = str(time_index).zfill(4)

        outpath = join(
            outdir,
            f"specriver_{ticid}_{linestr}_{tstr}_{lc_cadences}{s}.png"
        )

        fig.savefig(outpath, bbox_inches='tight', dpi=450)
        print(f"saved {outpath}")
