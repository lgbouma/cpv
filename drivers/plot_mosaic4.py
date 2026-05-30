"""
Four-panel horizontal mosaic: LP 12-502, DG CVn, TIC 300651846, TIC 262400835.
Light curves fetched via _get_lcpaths_given_ticid; periodsearch results cached in cachedir.
"""
import os
from os.path import join
from glob import glob

import numpy as np
import matplotlib.pyplot as plt

from astrobase.lcmath import time_bin_magseries

from aesthetic.plot import set_style, savefig

from complexrotators.paths import RESULTSDIR, LOCALDIR, LKCACHEDIR
from complexrotators.getters import _get_lcpaths_given_ticid
from complexrotators.lcprocessing import cpv_periodsearch, prepare_cpv_light_curve
from complexrotators.plotting import plot_phased_light_curve


def get_ylimguess(y):
    ylow = np.nanpercentile(y, 2.5)
    yhigh = np.nanpercentile(y, 97.5)
    ydiff = yhigh - ylow
    return [ylow - 0.35 * ydiff, yhigh + 0.35 * ydiff]


def subtract_secondary_sinusoid(time, flux, period_range_hr=(2.58, 2.62), n_periods=2000):
    """Grid-search for best secondary period in period_range_hr, subtract it.

    Fits y = A*sin(ωt) + B*cos(ωt) + C via linear least-squares on
    uncontaminated points (|flux-1| < 0.02), picks the period minimising
    the residual sum-of-squares, then removes only the oscillatory
    A*sin + B*cos component so the mean flux level is preserved.

    Returns (cleaned_flux, best_period_hr).
    """
    sel = np.abs(flux - 1) < 0.02
    t_fit, f_fit = time[sel], flux[sel]

    period_grid = np.linspace(period_range_hr[0] / 24, period_range_hr[1] / 24, n_periods)

    best_rss, best_p = np.inf, period_grid[0]
    for p in period_grid:
        omega = 2 * np.pi / p
        X = np.column_stack([np.sin(omega * t_fit),
                             np.cos(omega * t_fit),
                             np.ones_like(t_fit)])
        params, _, _, _ = np.linalg.lstsq(X, f_fit, rcond=None)
        rss = np.sum((f_fit - X @ params) ** 2)
        if rss < best_rss:
            best_rss = rss
            best_p = p

    omega = 2 * np.pi / best_p
    X_fit = np.column_stack([np.sin(omega * t_fit),
                             np.cos(omega * t_fit),
                             np.ones_like(t_fit)])
    a, b, _ = np.linalg.lstsq(X_fit, f_fit, rcond=None)[0]
    sinusoid = a * np.sin(omega * time) + b * np.cos(omega * time)
    best_p_hr = best_p * 24
    ampl_pct = np.sqrt(a**2 + b**2) * 1e2
    print(f'DG CVn secondary subtraction: best period = {best_p_hr:.4f} hr, amplitude = {ampl_pct:.3f}%')
    return flux - sinusoid, best_p_hr


# (ticid, sector, lcpipeline): LP 12-502, DG CVn, TIC 300651846, TIC 262400835
targets = [
    (300651846, 98, 'spoc2min'), # recent; good
    #(402980664, 85, 'spoc2min'), # only one dip!  slowly grows over sector
    (368129164, 50, 'spoc2min'), # good
    #(368129164, 77, 'spoc2min'), # good; & different!
    (402980664, 58, 'spoc2min'), # nice; complex
    #(300651846, 36, 'spoc2min'), # good
    (262400835, 32, 'spoc2min'), # good
    # (262400835, 98, 'qlp'), # good; recent (but QLP only)
]

PLOTDIR = join(RESULTSDIR, 'lc_mosaic')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

cachedir = join(LOCALDIR, 'cpv_finding', '2026_mosaic4')
if not os.path.exists(cachedir):
    os.makedirs(cachedir)

def plot_mosaic4():

    plt.close('all')
    set_style('science')

    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(1.2 * 4.67, 1.2 * 1.25),
                            constrained_layout=True)
    axs = axs.flatten()

    ticid_names = {
        '402980664': 'LP 12-502',
        '368129164': 'DG CVn',
        '300651846': 'TIC 300651846',
        '262400835': 'TIC 262400835',
    }

    for ix, ((ticid, sector, lcpipeline), ax) in enumerate(zip(targets, axs)):

        if lcpipeline in ['spoc2min', 'tess-spoc']:
            lcpaths = _get_lcpaths_given_ticid(str(ticid), lcpipeline)
        elif lcpipeline == 'qlp':
            lcpaths = glob(join(
                LKCACHEDIR.replace("TESS", "HLSP"),
                f'hlsp_qlp_tess_ffi_s{str(sector).zfill(4)}*{ticid}*llc',
                '*fits'
            ))

        sstr = str(sector).zfill(4)
        lcpath = [l for l in lcpaths if sstr in l][0]

        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector_out, starid) = prepare_cpv_light_curve(lcpath, cachedir)

        if str(ticid) == '368129164':
            y_flat, _ = subtract_secondary_sinusoid(x_obs, y_flat)

        d = cpv_periodsearch(
            x_obs, y_flat, starid, cachedir, t0='binmin', periodogram_method='pdm'
        )

        bd = time_bin_magseries(d['times'], d['fluxs'], binsize=200, minbinelems=1)
        ylim = get_ylimguess(1e2 * (bd['binnedmags'] - np.nanmean(bd['binnedmags'])))

        plot_phased_light_curve(
            d['times'], d['fluxs'], d['t0'], d['period'], None, ylim=ylim,
            xlim=[-0.6, 0.6], binsize_phase=1/300, BINMS=1.5, titlestr=None,
            showtext=None, showtitle=False, figsize=None, c0='darkgray',
            alpha0=0.3, c1='k', alpha1=1, phasewrap=True, plotnotscatter=False,
            fig=None, ax=ax, savethefigure=False, findpeaks_result=None,
            showxticklabels=False, titlefontsize='xx-small'
        )

        name = ticid_names.get(str(ticid), f'TIC {ticid}')
        titlestr = f'{name} (S{sector}, P={d["period"]*24:.2f} hr)'
        ax.set_title(titlestr, fontsize='xx-small', pad=2)

        ax.set_xticks([-0.5, 0, 0.5])

        ylow = int(np.ceil(ylim[0])) + 1
        yhigh = int(np.floor(ylim[1])) - 1
        if np.diff([np.abs(ylow), yhigh]) <= 2:
            ylowabs, yhighabs = np.abs(ylow), np.abs(yhigh)
            ylow = -np.min([ylowabs, yhighabs])
            yhigh = np.min([ylowabs, yhighabs])
            if ylow == yhigh == 0:
                ylow, yhigh = -1, 1
            if yhigh >= 10:
                ylow, yhigh = -9, 9

        ax.set_yticks([ylow, 0, yhigh])
        ax.set_yticklabels([ylow, 0, yhigh])

        id_yticks = {
            '402980664': [-2, 0, 2],
            '300651846': [-4, 0, 4],
            '262400835': [-4, 0, 4],
        }
        if str(ticid) in id_yticks:
            yticks = id_yticks[str(ticid)]
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks)

        if str(ticid) == '300651846':
            if sector == 98:
                ax.set_ylim([-6.25, 5.0])  # 25% bigger than [-5,4]
            else:
                ax.set_ylim([-5, 4])

        ax.tick_params(axis='both', which='major', labelsize='small')

    for ax in axs:
        ax.set_xticklabels(['-0.5', '0', '0.5'])

    fs = 'medium'
    axs[0].set_ylabel(r"$\Delta$ Flux [%]", fontsize=fs)
    fig.supxlabel(r"Phase, $\varphi$", fontsize=fs)

    outpath = join(PLOTDIR, 'lc_mosaic_mosaic4.png')
    savefig(fig, outpath, dpi=400)

if __name__ == "__main__":
    plot_mosaic4()
