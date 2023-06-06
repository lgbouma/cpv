"""
given a list of ticid's and their sectors, plot flux vs time for them, over
some number of cycles.
"""
#############
## IMPORTS ##
#############
import os, pickle
from os.path import join
from glob import glob
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from numpy import array as nparr

from aesthetic.plot import savefig, format_ax, set_style

from complexrotators.paths import LOCALDIR, SPOCDIR, RESULTSDIR
from complexrotators.getters import (
    _get_lcpaths_given_ticid, _get_local_lcpaths_given_ticid
)
from complexrotators.lcprocessing import (
    cpv_periodsearch, count_phased_local_minima, prepare_cpv_light_curve
)
from complexrotators.plotting import (
    plot_phased_light_curve, plot_dipcountercheck, plot_cpvvetter
)

from astrobase.lcmath import time_bin_magseries

from complexrotators import pipeline_utils as pu

############
# get data #
############

# define ticids and sectors

def plot_flux_vs_time(
    ticids_sectors = [
        (402980664,19),
        (94088626,10),
        (425933644,28),
    ],
    xmin=0,
    xmax=10
):

    #
    # get the light curves for all desired sectors and cadences
    #

    d = {}
    for ticid, sector in ticids_sectors:

        LOCAL_DEBUG = 1
        if LOCAL_DEBUG:
            lcpaths = _get_lcpaths_given_ticid(ticid)
        else:
            lcpaths = _get_local_lcpaths_given_ticid(ticid)

        lcpaths = [l for l in lcpaths if "s"+str(sector).zfill(4) in l]
        lcpath = lcpaths[0]

        # get the relevant light curve data
        cachedir = join(LOCALDIR, "flux_vs_time")
        if not os.path.exists(cachedir): os.mkdir(cachedir)

        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid) = prepare_cpv_light_curve(lcpath, cachedir)

        # get period, t0, and periodogram (PDM or LombScargle)
        _d = cpv_periodsearch(
            x_obs, y_flat, starid, cachedir, t0='binmin', periodogram_method='pdm'
        )

        # quantities we care about are d['times'], d['fluxs'], d['t0'], d['period']
        d[ticid] = _d


    #################
    # make the plot #
    #################
    outdir = join(RESULTSDIR, "flux_vs_time")
    if not os.path.exists(outdir): os.mkdir(outdir)

    N_rows = len(ticids_sectors)

    plt.close("all")
    set_style("clean")

    fig, axs = plt.subplots(nrows=N_rows, ncols=1, figsize=(1.2*3.5, 2*1.2*1.25),
                            constrained_layout=True, sharex=True)
    axs = axs.flatten()

    ny = lambda y: 1e2*(y-np.nanmean(y))

    for ix, (ticid, sector) in enumerate(ticids_sectors):

        P = d[ticid]['period']
        nx = lambda x: x-d[ticid]['t0']

        ax = axs[ix]

        # raw datta
        ax.scatter(nx(d[ticid]['times'])/P, ny(d[ticid]['fluxs']), c='lightgray', s=0.5,
                   linewidths=0, zorder=1, rasterized=True)

        # 1/100th of cycle bin.
        binsize_d = P/100
        binsize_sec = binsize_d*24*60*60
        print(f"TIC{ticid}, binsize {binsize_sec}s")
        bd = time_bin_magseries(d[ticid]['times'], d[ticid]['fluxs'],
                                binsize=binsize_sec, minbinelems=1)
        ax.plot(nx(bd['binnedtimes'])/P, ny(bd['binnedmags']), c='k', lw=0.5,
                zorder=2)

        ax.set_xlim([xmin, xmax])

        xsel = (nx(bd['binnedtimes'])/P > xmin) & (nx(bd['binnedtimes'])/P < xmax)
        ymin0 = min(ny(bd['binnedmags'])[xsel])
        ymax0 = max(ny(bd['binnedmags'])[xsel])
        yrange = ymax0 - ymin0
        ymin = ymin0 - 0.15*yrange
        ymax = ymax0 + 0.15*yrange
        ax.set_ylim([ymin, ymax])

        txt = f'$P$: {P*24:.1f} hr'
        props = dict(boxstyle='square', facecolor='white', alpha=0.7, pad=0.15,
                     linewidth=0)
        ax.text(0.99,0.03,txt,
                transform=ax.transAxes,
                ha='right',va='bottom', color='k', fontsize='small',
                bbox=props, zorder=3)


    fs = 'large'
    fig.text(-0.02,0.5, r"Flux [%]", va='center', rotation=90, fontsize=fs)
    axs[-1].set_xlabel(r"Time [cycle number]", fontsize=fs)

    # set naming options
    s = f'_xmin{str(xmin).zfill(4)}_xmax{str(xmax).zfill(4)}'

    outpath = join(outdir, f'flux_vs_time{s}.png')
    savefig(fig, outpath, dpi=400)

if __name__ == "__main__":
    plot_flux_vs_time(xmin=0, xmax=10)
    plot_flux_vs_time(xmin=10, xmax=20)
    plot_flux_vs_time(xmin=20, xmax=30)
