"""
Visualize the actual MagE spectra used by plot_mage_mosaic.py, one PNG per
timestamp, for a single night (default 2025/01/18).

For each of the night's exposures, make a three-panel figure of the Hα
echellette order (order index 3, ~6328-7363 Å) of TIC 300651846:

    A: full order, raw counts vs wavelength; Hα rest λ marked, and the
       ±DLAMBDA window used by the mosaic shaded.
    B: the ±DLAMBDA cutout, raw counts vs wavelength.
    C: the cutout exactly as plot_mage_mosaic.py uses it, i.e. Gaussian
       smoothed (σ=2 pix) and normalized to the smoothed continuum at
       NORMATVEL, plotted against Δv/v_eq.  The night's Gaussian-smoothed
       25th-percentile "average" template (the thing the specriver panels
       subtract) is shown in gray for reference.

Axis limits are held fixed across all timestamps so the PNGs flip through as a
movie.  A summary figure overplotting every spectrum, colored by phase, is also
written.

Conventions (period, T0, v_eq, dlambda, normalization, smoothing, order) are
inherited from drivers/plot_mage_mosaic.py and
complexrotators.getters.get_specriver_data.

env: cpv  (use /Users/luke/local/miniconda3/envs/cpv/bin/python)
"""

import os
from glob import glob
from os.path import join, basename

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u, constants as const
from astropy.time import Time
from scipy.ndimage import gaussian_filter1d

from aesthetic.plot import set_style, savefig

from complexrotators.paths import RESULTSDIR

# ── constants (matched to plot_mage_mosaic.py / get_specriver_data) ─────────
TICID = '300651846'
LINESTR = 'Hα'
LAMBDA0 = 6562.8                     # Hα rest wavelength [Å]
ORDER = 3                            # MagE order index holding Hα
START, END = 10, -10                 # order-edge pixel trim, as in getters.py
DLAMBDA = 15                         # cutout half-width [Å]
ADOPTED_VEQ = 89.65                  # km/s
NORMATVEL = 89.65 * 5.5              # km/s; normpointdict['300651846']
GAUSSIAN_FILTER_SIGMA = 2            # pix, applied to the cutout
PCTILE = 25                          # night-average template percentile
TEMPLATE_SMOOTH_SIGMA = 1            # pix, MaGE value used in the mosaic
PERIOD = 0.3439207894624283          # days
T0 = 2460693.243837336 - 2457000     # BTJD

DATESTR = '20250118'
DATE_LABEL = '2025/01/18'

SPECDIR = os.path.expanduser(
    '~/Dropbox/proj/cpv/data/spectra/MAGE_REDUX/Final-Products'
)
TIMECSV = join(SPECDIR, 'tic3006_mage_bjd_times.csv')

# glob patterns per night, copied from get_specriver_data
GLOBDICT = {
    '20250118': 'tic3006_spec0*.fits',
    '20250119': 'tic3006_spec1*.fits',
    '20250120': 'tic3006_spec2*.fits',
    '20250121': 'tic3006_spec3*.fits',
}
# 20250119 drops its first spectrum in get_specriver_data
DROPFIRST = {'20250119'}


def _get_vel(wav, wav0=LAMBDA0):
    return (const.c * ((wav - wav0) / wav0)).to(u.km / u.s).value


def load_night(datestr=DATESTR):
    """Load the Hα order for every exposure of one night.

    Returns a dict with the full order (wav_full/flx_full), the ±DLAMBDA
    cutout (wav_cut/flx_cut), the mosaic's normalized+smoothed cutout
    (flx_norm) and its Δv/v_eq abscissa (vel_cut), plus times/phases.
    """
    from cdips_followup.spectools import read_mage

    specpaths = np.sort(glob(join(SPECDIR, GLOBDICT[datestr])))
    if datestr in DROPFIRST:
        specpaths = specpaths[1:]
    assert len(specpaths) > 0

    timedf = pd.read_csv(TIMECSV)
    sel = timedf.observation_id.astype(str).str.contains(datestr)
    seltimedf = timedf[sel]
    spectimes = np.array(seltimedf.bjd_tdb_midtime) - 2457000
    dateobs = np.array(seltimedf.date_obs)
    exptimes = np.array(seltimedf.exptime_s)
    # same trailing-alignment guard as plot_mage_mosaic._make_specriver_data
    if len(spectimes) > len(specpaths):
        n_extra = len(spectimes) - len(specpaths)
        spectimes, dateobs, exptimes = (
            spectimes[n_extra:], dateobs[n_extra:], exptimes[n_extra:]
        )
    assert len(spectimes) == len(specpaths), (
        f'{len(spectimes)} times vs {len(specpaths)} spectra'
    )

    wav_full, flx_full, wav_cut, flx_cut, vel_cut, flx_norm, norm_flxs = (
        [], [], [], [], [], [], []
    )
    for ix, specpath in enumerate(specpaths):
        print(f'{ix} {basename(specpath)}')
        flx_2d, wav_2d = read_mage(specpath)
        flx, wav = flx_2d[ORDER, START:END], wav_2d[ORDER, START:END]
        wav_full.append(wav)
        flx_full.append(flx)

        m = (wav > LAMBDA0 - DLAMBDA) & (wav < LAMBDA0 + DLAMBDA)
        _w, _f = wav[m], flx[m]
        vels = _get_vel(_w)
        # normalize at NORMATVEL off the σ=5 smoothed cutout, as in getters.py
        norm_flx = gaussian_filter1d(_f, sigma=5)[np.argmin(abs(vels - NORMATVEL))]
        wav_cut.append(_w)
        flx_cut.append(_f)
        vel_cut.append(vels / ADOPTED_VEQ)
        flx_norm.append(
            gaussian_filter1d(_f / norm_flx, sigma=GAUSSIAN_FILTER_SIGMA)
        )
        norm_flxs.append(norm_flx)

    # phases in the mosaic's per-night convention (first exposure near φ=0)
    phases = (spectimes - T0) / PERIOD
    phases -= np.min(np.ceil(phases))

    # night-average template: Gaussian-smoothed PCTILE-th percentile spectrum.
    # Cutouts can differ by one pixel across the night (wavelength solution
    # drift); trim to the shortest, as get_specriver_data does.
    nmin = min(len(f) for f in flx_norm)
    stack = np.array([f[:nmin] for f in flx_norm]).T   # (n_wav, n_spec)
    template = gaussian_filter1d(
        np.nanpercentile(stack, PCTILE, axis=1), sigma=TEMPLATE_SMOOTH_SIGMA
    )

    return {
        'specpaths': specpaths,
        'spectimes': spectimes,
        'dateobs': dateobs,
        'exptimes': exptimes,
        'phases': phases,
        'wav_full': wav_full,
        'flx_full': flx_full,
        'wav_cut': wav_cut,
        'flx_cut': flx_cut,
        'vel_cut': vel_cut,
        'flx_norm': flx_norm,
        'norm_flxs': np.array(norm_flxs),
        'template': template,
        'template_vel': vel_cut[0][:nmin],
    }


def _shared_limits(d):
    """Fixed axis limits across the night, so frames are comparable."""
    all_full = np.concatenate(d['flx_full'])
    all_cut = np.concatenate(d['flx_cut'])
    all_norm = np.concatenate(d['flx_norm'])
    # top of the full-order axis must clear the strongest Hα peak of the night
    # (all_cut.max()) while staying robust to hot pixels elsewhere in the order
    hi_full = max(np.nanpercentile(all_full, 99.9), np.nanmax(all_cut))
    hi_cut = np.nanpercentile(all_cut, 99.9)
    return {
        'ylim_full': (0, 1.08 * hi_full),
        'ylim_cut': (0, 1.12 * hi_cut),
        'ylim_norm': (
            np.nanpercentile(all_norm, 0.2) - 0.1,
            np.nanpercentile(all_norm, 99.8) + 0.1,
        ),
        'xlim_full': (
            min(w.min() for w in d['wav_full']),
            max(w.max() for w in d['wav_full']),
        ),
        'xlim_cut': (LAMBDA0 - DLAMBDA, LAMBDA0 + DLAMBDA),
        'xlim_vel': (
            _get_vel(LAMBDA0 - DLAMBDA) / ADOPTED_VEQ,
            _get_vel(LAMBDA0 + DLAMBDA) / ADOPTED_VEQ,
        ),
    }


def _plot_one(ix, d, lim, outdir, datestr=DATESTR, datelabel=DATE_LABEL):

    plt.close('all')
    set_style('science')

    fig = plt.figure(figsize=(7, 4.2))
    axd = fig.subplot_mosaic(
        """
        AA
        BC
        """,
        gridspec_kw={'height_ratios': [1, 1], 'hspace': 0.45, 'wspace': 0.28},
    )

    wav, flx = d['wav_full'][ix], d['flx_full'][ix]
    t = d['spectimes'][ix]
    phase = d['phases'][ix]
    ut = str(d['dateobs'][ix]).replace('T', ' ')

    # ── A: full Hα order ───────────────────────────────────────────────────
    ax = axd['A']
    ax.plot(wav, flx, c='k', lw=0.3, rasterized=True)
    ax.axvline(LAMBDA0, c='crimson', lw=0.5, ls='--', alpha=0.7, zorder=1)
    ax.fill_between(
        [LAMBDA0 - DLAMBDA, LAMBDA0 + DLAMBDA], *lim['ylim_full'],
        color='crimson', alpha=0.12, linewidth=0, zorder=0,
    )
    ax.set_xlim(lim['xlim_full'])
    ax.set_ylim(lim['ylim_full'])
    ax.set_xlabel(r'$\lambda$ [$\AA$]', fontsize='x-small')
    ax.set_ylabel(r'$f_\lambda$ [counts]', fontsize='x-small')
    ax.tick_params(labelsize='x-small')
    ax.text(
        0.995, 0.95, f'order {ORDER}', ha='right', va='top',
        transform=ax.transAxes, fontsize='xx-small',
    )
    txt = (
        f'{datelabel}  |  {ix:02d}/{len(d["specpaths"]) - 1}  |  '
        f'{basename(d["specpaths"][ix]).replace("ofso_multi.fits", "")}\n'
        f'BTJD {t:.5f}  |  UT {ut}  |  '
        f'{d["exptimes"][ix]:.0f} s  |  '+ r'$\varphi$ = ' + f'{phase:.3f}'
    )
    ax.set_title(txt, fontsize='xx-small', pad=4)

    # ── B: raw cutout ──────────────────────────────────────────────────────
    ax = axd['B']
    ax.plot(d['wav_cut'][ix], d['flx_cut'][ix], c='k', lw=0.6)
    ax.axvline(LAMBDA0, c='crimson', lw=0.5, ls='--', alpha=0.7)
    ax.set_xlim(lim['xlim_cut'])
    ax.set_ylim(lim['ylim_cut'])
    ax.set_xlabel(r'$\lambda$ [$\AA$]', fontsize='x-small')
    ax.set_ylabel(r'$f_\lambda$ [counts]', fontsize='x-small')
    ax.tick_params(labelsize='x-small')
    ax.text(
        0.04, 0.95, f'{LINESTR} ' + r'$\pm$' + f'{DLAMBDA}' + r'$\,\AA$',
        ha='left', va='top', transform=ax.transAxes, fontsize='xx-small',
        bbox={'facecolor': 'white', 'alpha': 0.5, 'edgecolor': 'none',
              'pad': 0.4},
    )

    # ── C: normalized + smoothed cutout, as used by the mosaic ─────────────
    ax = axd['C']
    ax.plot(d['template_vel'], d['template'], c='darkgray', lw=0.8, zorder=1,
            label=f'night {PCTILE}th pctile')
    ax.plot(d['vel_cut'][ix], d['flx_norm'][ix], c='k', lw=0.8, zorder=2,
            label='this spectrum')
    ax.axvline(0, c='crimson', lw=0.5, ls='--', alpha=0.7, zorder=0)
    for v in (-1, 1):
        ax.axvline(v, c='darkgray', lw=0.4, ls=':', alpha=0.8, zorder=0)
    ax.set_xlim(lim['xlim_vel'])
    ax.set_ylim(lim['ylim_norm'])
    ax.set_xlabel(r'$\Delta v/v_\mathrm{eq}$', fontsize='x-small')
    ax.set_ylabel(r'$f_\lambda$ / $f_\lambda(v_\mathrm{norm})$',
                  fontsize='x-small')
    ax.tick_params(labelsize='x-small')
    ax.legend(fontsize='xx-small', loc='upper right', framealpha=0.8)

    fname = basename(d['specpaths'][ix]).replace('ofso_multi.fits', '')
    outpath = join(outdir, f'mage_{datestr}_order{ORDER}_{ix:02d}_{fname}.png')
    savefig(fig, outpath, writepdf=False, dpi=200)


def _plot_summary(d, lim, outdir, datestr=DATESTR, datelabel=DATE_LABEL):
    """All spectra of the night overplotted, colored by phase."""
    plt.close('all')
    set_style('science')

    fig, axs = plt.subplots(nrows=2, figsize=(7, 5))
    cmap = mpl.cm.viridis
    norm = mpl.colors.Normalize(
        vmin=d['phases'].min(), vmax=d['phases'].max()
    )

    for ix in range(len(d['specpaths'])):
        c = cmap(norm(d['phases'][ix]))
        axs[0].plot(d['wav_full'][ix], d['flx_full'][ix], c=c, lw=0.2,
                    alpha=0.7, rasterized=True)
        axs[1].plot(d['vel_cut'][ix], d['flx_norm'][ix], c=c, lw=0.5,
                    alpha=0.8)

    axs[0].axvline(LAMBDA0, c='crimson', lw=0.5, ls='--', alpha=0.7)
    axs[0].set_xlim(lim['xlim_full'])
    axs[0].set_ylim(lim['ylim_full'])
    axs[0].set_xlabel(r'$\lambda$ [$\AA$]', fontsize='x-small')
    axs[0].set_ylabel(r'$f_\lambda$ [counts]', fontsize='x-small')
    axs[0].set_title(
        f'{datelabel}: MagE order {ORDER}, N={len(d["specpaths"])} spectra',
        fontsize='x-small', pad=4,
    )

    axs[1].plot(d['template_vel'], d['template'], c='k', lw=1.0, zorder=5)
    axs[1].axvline(0, c='crimson', lw=0.5, ls='--', alpha=0.7)
    axs[1].set_xlim(lim['xlim_vel'])
    axs[1].set_ylim(lim['ylim_norm'])
    axs[1].set_xlabel(r'$\Delta v/v_\mathrm{eq}$', fontsize='x-small')
    axs[1].set_ylabel(r'$f_\lambda$ / $f_\lambda(v_\mathrm{norm})$',
                      fontsize='x-small')

    for ax in axs:
        ax.tick_params(labelsize='x-small')

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cax = axs[1].inset_axes([1.02, 0.0, 0.03, 1.0])
    cb = fig.colorbar(sm, cax=cax, orientation='vertical')
    cb.set_label(r'Phase, $\varphi$', rotation=90, fontsize='xx-small')
    cb.ax.tick_params(labelsize='xx-small')

    outpath = join(outdir, f'mage_{datestr}_order{ORDER}_summary.png')
    savefig(fig, outpath, writepdf=False, dpi=250)


def plot_mage_halpha_order_viz(datestr=DATESTR, datelabel=DATE_LABEL):

    outdir = join(RESULTSDIR, 'mage_viz', f'{datestr}_halpha_order')
    os.makedirs(outdir, exist_ok=True)

    d = load_night(datestr)
    lim = _shared_limits(d)

    for ix in range(len(d['specpaths'])):
        _plot_one(ix, d, lim, outdir, datestr=datestr, datelabel=datelabel)

    _plot_summary(d, lim, outdir, datestr=datestr, datelabel=datelabel)

    print(f'\nWrote {len(d["specpaths"]) + 1} PNGs to {outdir}')


if __name__ == '__main__':
    plot_mage_halpha_order_viz()
