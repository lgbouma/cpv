import os
from os.path import join
from glob import glob

import numpy as np, matplotlib.pyplot as plt, pandas as pd

from aesthetic.plot import set_style, savefig

from complexrotators.getters import (
    _get_lcpaths_fromlightkurve_given_ticid, get_qlp_lcpaths
)
from complexrotators.lcprocessing import (
    cpv_periodsearch, prepare_cpv_light_curve
)
from complexrotators.paths import RESULTSDIR, DATADIR, LKCACHEDIR
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




def make_plot(ticid, sector=None, showtitles=0, showphase=1,
              lcpipeline='spoc2min', style='clean', xlim=None,
              do_quality_trim=1, period=None, edge_trim=0, bincadence=None):

    # get data
    if lcpipeline in ['spoc2min', 'tess-spoc']:
        lcpaths = _get_lcpaths_fromlightkurve_given_ticid(ticid, lcpipeline)
    elif lcpipeline == 'qlp':
        lcpaths = glob(join(
            LKCACHEDIR.replace("TESS","HLSP"),
            f'hlsp_qlp_tess_ffi_s{str(sector).zfill(4)}*{ticid}*llc',
            '*fits'
        ))
        assert len(lcpaths) > 0

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
         lcpath, cachedir, lcpipeline=lcpipeline,
         do_quality_trim=do_quality_trim
     )

    # get period, t0, and periodogram (PDM or LombScargle)
    d = cpv_periodsearch(
        x_obs, y_flat, starid, cachedir, t0='binmin', periodogram_method='pdm'
    )
    if period is not None:
        d['period'] = period

    bd = time_bin_magseries(x_obs, y_flat, binsize=1200, minbinelems=1)
    #if not edge_trim:
    #    bd = time_bin_magseries(d['times'], d['fluxs'], binsize=1200, minbinelems=1)
    #else:
    #    _x = d['times']
    #    _sel = (_x - np.nanmin(_x) > 1) & (_x - np.nanmin(_x) < 25)
    #    bd = time_bin_magseries(d['times'][_sel], d['fluxs'][_sel], binsize=1200, minbinelems=1)

    #ylim = get_ylimguess(1e2*(bd['binnedmags']-np.nanmean(bd['binnedmags'])))
    ylim = get_ylimguess(y_flat)

    if showtitles:
        titlestr = f'{ticid}, s{sector}, {d["period"]*24:.1f}h'
    else:
        titlestr = None

    binsize_phase = 1/150
    BINMS=1.5
    alpha0=0.3

    # make plot
    plt.close('all')
    set_style(style)
    if showphase:
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(0.6666*1.2*4.5, 1.2*1.25),
                                constrained_layout=True)
        axs = axs.flatten()
    else:
        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(0.6666*1.2*4.5, 1.2*1.25),
                                constrained_layout=True)
        axs = [axs]

    c = 'k' if 'wob' not in style else 'white'
    c2 = 'lightgray'


    if bincadence is None:
        x_offset = np.nanmin(x_obs)
        if ticid == '314847177':
            x_offset = np.nanmin(x_obs) + 15
        axs[0].scatter(x_obs-x_offset, y_flat, c=c, s=0.5, linewidths=0, zorder=10)
        _, _groups = find_lc_timegroups(x_obs, mingap=0.5/24)
        for _g in _groups:
            axs[0].plot(x_obs[_g]-x_offset, y_flat[_g], c=c2, zorder=8, lw=0.2, alpha=0.25)
    else:
        _bd = time_bin_magseries(x_obs, y_flat, binsize=bincadence, minbinelems=1)
        _x, _y = _bd['binnedtimes'], _bd['binnedmags']
        _x_offset = np.nanmin(_x)
        axs[0].scatter(_x-_x_offset, _y, c=c, s=0.5, linewidths=0, zorder=10)
        _, _groups = find_lc_timegroups(_x, mingap=0.5/24)
        for _g in _groups:
            axs[0].plot(_x[_g]-_x_offset, _y[_g], c=c2, zorder=8, lw=0.2, alpha=0.25)


    axs[0].set_xlabel('Time [days]')
    axs[0].set_ylabel('Relative flux')
    if showphase:
        axs[1].set_ylabel('Relative flux')

    axs[0].set_xlim([-0.2,10.2])
    if isinstance(xlim, (list, tuple)):
        axs[0].set_xlim(xlim)
    if ticid == '220599904':
        axs[0].set_xlim([14, 24])

    if showphase:
        c1 = 'k' if 'wob' not in style else 'white'
        plot_phased_light_curve(
            d['times'], d['fluxs'], d['t0'], d['period'], None, ylim=ylim,
            xlim=[-0.6,0.6], binsize_phase=binsize_phase, BINMS=BINMS, titlestr=titlestr,
            showtext=False, showtitle=False, figsize=None, c0='darkgray',
            alpha0=alpha0, c1=c1, alpha1=1, phasewrap=True, plotnotscatter=False,
            fig=None, ax=axs[1], savethefigure=False, findpeaks_result=None,
            showxticklabels=[-0.5,0,0.5], titlefontsize=0, normfunc=0
        )
        txt = f'$P$={d["period"]*24:.1f}$\,$hr'
        axs[1].text(0.97,0.97,txt, transform=axs[1].transAxes,
                    ha='right',va='top')

        axs[1].update({'xlabel': 'Phase, φ'})

    for ax in axs:
        ax.set_ylim(ylim)

    outdir = join(RESULTSDIR, "normal_stars")
    if not os.path.exists(outdir): os.mkdir(outdir)

    s = '' if 'wob' not in style else '_wob'
    p = '' if showphase else '_nophase'
    outpath = join(outdir, f"TIC{ticid}_s{str(sector).zfill(4)}{s}{p}.png")
    savefig(fig, outpath, dpi=400)

if __name__ == "__main__":

    styles = ['clean_wob', 'clean']
    for style in styles:

        # AB Dor, great rotator
        make_plot('149248196', sector=7, style=style)
        make_plot('149248196', sector=7, style=style, showphase=0,
                  xlim=[-0.2, 6.2])
        assert 0

        # OO Peg, great EB
        make_plot('314847177', sector=55, lcpipeline='tess-spoc', style=style,
                 xlim=[-0.2, 28.2], do_quality_trim=0, period=2.98465593,
                  edge_trim=1)
        make_plot('314847177', sector=55, lcpipeline='tess-spoc', style=style,
                  showphase=0, xlim=[-0.2, 12.2], do_quality_trim=0,
                  period=2.98465593, edge_trim=1)

        # LP 12-502 cpv
        #for sector in [18, 19, 25, 26, 53, 58, 73]:
        for sector in [58]:
            make_plot("402980664", sector=sector, style=style)
            make_plot("402980664", sector=sector, style=style, showphase=0,
                      bincadence=600)
        # disk
        make_plot('57528302', sector=36, showphase=0, style=style)

        # nice rotator (secretly cpv sometimes)
        make_plot('177309964', sector=67, style=style)
        make_plot('177309964', sector=67, style=style, showphase=0)

        ## nice eb
        #make_plot('281498280', sector=1, style=style)
        #make_plot('281498280', sector=1, style=style, showphase=0)

        assert 0
        # qlp cpv
        make_plot("260268310", sector=31, lcpipeline='qlp', showphase=1)

        # qlp cpv
        make_plot("35858638", sector=31, lcpipeline='qlp', showphase=1)

        # qlp cpv
        make_plot("220599904", sector=31, lcpipeline='qlp', showphase=0)

        # AU Mic
        make_plot('441420236', sector=1)

        # HIP 67522
        make_plot("166527623")

    # # # shape-changing rotator
    # make_plot("93839949")
