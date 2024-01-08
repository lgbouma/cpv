import os
from os.path import join

import numpy as np, matplotlib.pyplot as plt, pandas as pd

from aesthetic.plot import set_style, savefig

from complexrotators.getters import (
    _get_lcpaths_fromlightkurve_given_ticid, get_qlp_lcpaths
)
from complexrotators.lcprocessing import (
    cpv_periodsearch, prepare_cpv_light_curve
)
from complexrotators.paths import RESULTSDIR, DATADIR
from complexrotators.plotting import plot_phased_light_curve

from astrobase.lcmath import (
    phase_magseries, phase_bin_magseries, sigclip_magseries,
    find_lc_timegroups, phase_magseries_with_errs, time_bin_magseries
)

def get_ylimguess(y):
    ylow = np.nanpercentile(y, 2.5)
    yhigh = np.nanpercentile(y, 97.5)
    ydiff = (yhigh-ylow)
    ymin = ylow - 0.35*ydiff
    ymax = yhigh + 0.35*ydiff
    return [ymin,ymax]




def make_plot(ticid, sector=None, showtitles=0, showphase=1, lcpipeline='spoc2min'):

    # get data
    if lcpipeline == 'spoc2min':
        lcpaths = _get_lcpaths_fromlightkurve_given_ticid(ticid)
    elif lcpipeline == 'qlp':
        # FIXME temp hack...
        lcpaths = [join(DATADIR, 'photometry', 'tess',
                        'hlsp_qlp_tess_ffi_s0031-0000000220599904_tess_v01_llc.fits')]
        #lcpaths = get_qlp_lcpaths(ticid)

    if isinstance(sector,int):
        sstr = str(sector).zfill(4)
        lcpath = [l for l in lcpaths if sstr in l][0]
    else:
        lcpath = lcpaths[0]

    cachedir = '/Users/luke/local/complexrotators/cpv_finding/normal_stars'
    if not os.path.exists(cachedir): os.mkdir(cachedir)

    # get the relevant light curve data
    (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
     sector, starid) = prepare_cpv_light_curve(
         lcpath, cachedir, lcpipeline=lcpipeline
     )

    # get period, t0, and periodogram (PDM or LombScargle)
    d = cpv_periodsearch(
        x_obs, y_flat, starid, cachedir, t0='binmin', periodogram_method='pdm'
    )

    bd = time_bin_magseries(d['times'], d['fluxs'], binsize=1200, minbinelems=1)
    #ylim = get_ylimguess(1e2*(bd['binnedmags']-np.nanmean(bd['binnedmags'])))
    ylim = get_ylimguess(y_flat)

    if showtitles:
        titlestr = f'{ticid}, s{sector}, {d["period"]*24:.1f}h'
    else:
        titlestr = None

    binsize_phase = 1/300
    BINMS=1.5
    alpha0=0.3

    # make plot
    plt.close('all')
    set_style('clean')
    if showphase:
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(0.6666*1.2*3.5, 1.2*1.25),
                                constrained_layout=True, sharey=True)
        axs = axs.flatten()
    else:
        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(0.6666*1.2*3.5, 1.2*1.25),
                                constrained_layout=True, sharey=True)
        axs = [axs]

    axs[0].scatter(x_obs-np.nanmin(x_obs)-1, y_flat, c='k', s=0.35, linewidths=0, zorder=10)
    axs[0].set_xlabel('Time [days]')
    axs[0].set_ylabel('Relative flux')

    if not ticid == '220599904':
        axs[0].set_xlim([-1,11])
    else:
        axs[0].set_xlim([14, 24])

    if showphase:
        plot_phased_light_curve(
            d['times'], d['fluxs'], d['t0'], d['period'], None, ylim=ylim,
            xlim=[-0.6,0.6], binsize_phase=binsize_phase, BINMS=BINMS, titlestr=titlestr,
            showtext=False, showtitle=False, figsize=None, c0='darkgray',
            alpha0=alpha0, c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
            fig=None, ax=axs[1], savethefigure=False, findpeaks_result=None,
            showxticklabels=[-0.5,0,0.5], titlefontsize=0, normfunc=0
        )
        txt = f'P={d["period"]:.3f}d'
        axs[1].text(0.97,0.97,txt, transform=axs[1].transAxes, ha='right',va='top',
                    color='k', fontsize='x-small')

        axs[1].update({'xlabel': 'Phase, φ'})

    outdir = join(RESULTSDIR, "normal_stars")
    if not os.path.exists(outdir): os.mkdir(outdir)

    outpath = join(outdir, f"TIC{ticid}_s{str(sector).zfill(4)}.png")
    savefig(fig, outpath, dpi=400)

if __name__ == "__main__":

    # qlp cpv
    make_plot("220599904", sector=31, lcpipeline='qlp', showphase=0)

    # AU Mic
    make_plot('441420236', sector=1)

    # disk
    make_plot('57528302', sector=36, showphase=0)

    # nice rotator (secretly cpv sometimes)
    make_plot('177309964', sector=67)

    # nice eb
    make_plot('281498280', sector=1)

    # HIP 67522
    make_plot("166527623")


    # # # shape-changing rotator
    # make_plot("93839949")


