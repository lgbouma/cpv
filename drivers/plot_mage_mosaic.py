"""
MaGE mosaic: TIC 300651846 four-night spectroscopic overview.

Layout:
    AAAB
    AAAB
    CDEF
    GHIJ
    GHIJ

A: normalized flux vs time, sector 88, xlim [3694.25, 3698.5];
   green bars for spectral coverage per night.
B: normalized flux vs phase, sector 88, full sector duration.
C/D/E/F: flux vs phase + Hα linecore twinx (= sixpanel panel D),
         for dates 20250118 / 20250119 / 20250120 / 20250121.
G/H/I/J: specriver pcolor removeavg (= sixpanel panel E),
         for dates 20250118 / 20250119 / 20250120 / 20250121.

env: cpv  (use /Users/luke/local/miniconda3/envs/cpv/bin/python)
"""

import os
from os.path import join
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.ndimage import gaussian_filter1d

from astrobase.lcmath import phase_magseries, phase_bin_magseries

from aesthetic.plot import set_style, savefig

from complexrotators.paths import RESULTSDIR, LOCALDIR
from complexrotators.getters import get_specriver_data
from complexrotators.plotting import plot_phased_light_curve, prepare_local_lc
from complexrotators.lcprocessing import cpv_periodsearch

# ── constants ──────────────────────────────────────────────────────────────
TICID = '300651846'
LINESTR = 'Hα'
DLAMBDA = 15
PERIOD = 0.3439207894624283          # days
T0 = 2460693.243837336 - 2457000     # TJD (BTJD), same as sixpanel driver
YLIM = [-8, 5]                       # flux vs phase ylim [%]
ADOPTED_VEQ = 89.65                  # km/s
PCTILE = 25                          # for removeavg gaussian smoothed percentile
BINSIZE_PHASE = 0.005
VMIN_RIVER, VMAX_RIVER = -0.2, 0.6  # specriver colormap range

DATES = ['20250118', '20250119', '20250120', '20250121']
DATE_LABELS = ['2025/01/18', '2025/01/19', '2025/01/20', '2025/01/21']

LCDIR = join(os.path.expanduser('~'), 'Dropbox/proj/cpv/data/photometry/tess')
LCPATH = join(LCDIR, 'tess2025014115807-s0088-0000000300651846-0285-s_lc.fits')
CACHEDIR = join(LOCALDIR, 'cpv_finding', 'spoc2min_debug')

PLOTDIR = join(RESULTSDIR, 'mage_mosaic')
os.makedirs(PLOTDIR, exist_ok=True)


def _make_specriver_data(datestr):
    """Load MaGE spectra for one night and return processed arrays.

    Returns dict with keys:
        spectimes   : BJD-2457000 midtimes (n_spec,)
        specphases  : phases centred near 0 for that night
        xvals       : list of Δv/v_eq arrays, one per spectrum
        orig_yvals  : list of normalized flux arrays before avg subtraction
        flux_arr_sub: (n_wav, n_spec) flux minus smoothed average
        norm_flxs   : continuum normalization values (n_spec,)
        linecore_sums: sum of line-core flux for each spectrum
    """
    specpaths, spectimes, xvals, yvals, yvalsnonorm, norm_flxs = get_specriver_data(
        TICID, LINESTR, dlambda=DLAMBDA, usespectype='mage', datestr=datestr
    )

    # normalize x to Δv/v_eq
    xvals = [np.asarray(xv) / ADOPTED_VEQ for xv in xvals]

    # spectimes may have one extra leading entry vs yvals (e.g. 20250119 has an
    # excluded first spectrum in specpaths but its time is still in the CSV).
    spectimes = np.asarray(spectimes)
    n_spec = len(yvals)
    if len(spectimes) > n_spec:
        spectimes = spectimes[len(spectimes) - n_spec:]

    # phases centred so that the first observation starts near 0
    xval_phase = (spectimes - T0) / PERIOD
    xval_phase -= np.min(np.ceil(xval_phase))

    # build (n_wav × n_spec) flux array
    n_wav = len(xvals[0])
    flux_arr = np.zeros((n_wav, n_spec))
    for ix, yval in enumerate(yvals):
        flux_arr[:, ix] = yval

    # linecore sums before avg subtraction
    linecore_sums = np.array([
        np.sum(yval[np.abs(xval) <= 1])
        for yval, xval in zip(yvals, xvals)
    ])

    # removeavg: subtract gaussian-smoothed 25th-percentile spectrum
    pctflux = np.nanpercentile(flux_arr, PCTILE, axis=1)
    smoothmeanflux = gaussian_filter1d(pctflux, sigma=1)  # sigma=1 for MaGE

    orig_yvals = deepcopy(yvals)
    for yval in yvals:
        yval -= smoothmeanflux
    flux_arr_sub = flux_arr - smoothmeanflux[:, None]

    return {
        'spectimes': np.asarray(spectimes),
        'specphases': xval_phase,
        'xvals': xvals,
        'orig_yvals': orig_yvals,
        'flux_arr_sub': flux_arr_sub,
        'norm_flxs': norm_flxs,
        'linecore_sums': linecore_sums,
    }


def _plot_phase_linecore(ax, lc_times, lc_fluxs, sd, datestr):
    """Draw flux vs phase + Hα linecore twinx (sixpanel panel D) on ax.

    Replicates the axd['D'] + ax2 linecore logic from
    plot_movie_sixpanel_specriver, style='science' (black on white).

    Returns the specriver phase limits (xmin, xmax) for use in panel G/H/I/J.
    """
    # ── linecore twinx (drawn first so TESS LC renders on top) ─────────────
    linecore_sums = sd['linecore_sums']
    norm_flxs = sd['norm_flxs']

    linecore_rel = linecore_sums / np.nanmedian(linecore_sums)
    med_normflx = np.nanmedian(norm_flxs)
    linecore_err = (
        np.sqrt(linecore_sums * med_normflx) /
        (linecore_sums * med_normflx)
    )
    linecore_rel_pct = 100.0 * (linecore_rel - np.nanmedian(linecore_rel))
    linecorr_err_pct = 100.0 * linecore_err
    FUDGE = 2

    specphases = sd['specphases']
    xerr = np.nanmedian(np.diff(sd['spectimes']) / PERIOD) / 2

    ax2 = ax.twinx()
    ax2.errorbar(
        specphases, linecore_rel_pct,
        yerr=FUDGE * linecorr_err_pct, xerr=xerr,
        lw=0.8, ls=':', marker='.', c='forestgreen', markersize=2,
        alpha=0.5, zorder=1,
    )
    ax2.set_ylabel(r"$\Delta$ $f$$_{\mathrm{H\alpha\ core}}$ [%]",
                   fontsize='x-small', color='forestgreen')
    ax2.tick_params(axis='y', labelcolor='forestgreen', labelsize='x-small')

    # cap outlier-driven scale for 20250119
    if datestr == '20250119':
        ax2.set_ylim(-25, 25)

    # bring TESS LC axis in front of the twinx linecore axis
    ax2.set_zorder(1)
    ax.set_zorder(2)
    ax.patch.set_visible(False)   # keep ax2 background visible

    # ── phased LC ──────────────────────────────────────────────────────────
    plot_phased_light_curve(
        lc_times, lc_fluxs, T0, PERIOD, None,
        fig=ax.get_figure(), ax=ax,
        binsize_phase=BINSIZE_PHASE,
        xlim=None,
        ylim=YLIM,
        showtext=None,
        showtitle=False,
        savethefigure=False,
        c0='darkgray', alpha0=0.15,
        c1='k', alpha1=1,
        phasewrap=False,
        longwrap=True,
        BINMS=1.5,
    )
    ax.set_ylim(YLIM)

    xmin = specphases.min() - 0.05
    xmax = specphases.max() + 0.05
    ax.set_xlim(xmin, xmax)

    ax.set_ylabel(r"$\Delta$ Flux [%]", fontsize='x-small')
    ax.set_xlabel(r"Phase, $\varphi$", fontsize='x-small')
    ax.tick_params(labelsize='x-small')

    return xmin, xmax, ax2


def _plot_specriver(ax, fig, sd):
    """Draw specriver pcolor removeavg (sixpanel panel E) on ax.

    Returns the pcolor mappable for colorbar use, and (xmin, xmax).
    """
    xvals = sd['xvals']     # list of Δv/v_eq arrays
    specphases = sd['specphases']
    flux_arr_sub = sd['flux_arr_sub']

    # In vertphase orientation:
    #   x-axis = phase of observations
    #   y-axis = Δv/v_eq
    #   color  = flux_arr_sub (n_wav × n_spec)
    y_vel = xvals[0]        # Δv/v_eq (assumed same for all spectra)
    x_phase = specphases    # phases (n_spec,)

    norm = colors.Normalize(vmin=VMIN_RIVER, vmax=VMAX_RIVER)
    c = ax.pcolor(
        x_phase, y_vel, flux_arr_sub,
        cmap='Spectral_r', norm=norm, shading='auto', rasterized=True,
    )

    ax.set_ylabel(r"$\Delta v/v_\mathrm{eq}$", fontsize='x-small')
    ax.set_xlabel(r"Phase, $\varphi$", fontsize='x-small')
    ax.tick_params(labelsize='x-small')

    # label "Hα - Avg."
    ax.text(
        0.04, 0.95, f'{LINESTR} - Avg.',
        ha='left', va='top', transform=ax.transAxes,
        fontsize='xx-small',
        bbox={"facecolor": "white", "alpha": 0.5, "edgecolor": "none", "pad": 0.4},
    )

    xmin, xmax = ax.get_xlim()
    return c, xmin, xmax


def plot_mage_mosaic():

    # ── load LC ────────────────────────────────────────────────────────────
    (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend,
     cadence_sec, sector, starid) = prepare_local_lc(TICID, LCPATH, CACHEDIR)

    d_ps = cpv_periodsearch(x_obs, y_flat, starid, CACHEDIR, t0='binmin')
    lc_times = d_ps['times']
    lc_fluxs = d_ps['fluxs']

    # ── load spectral data ─────────────────────────────────────────────────
    print('Loading spectral data...')
    spec = {d: _make_specriver_data(d) for d in DATES}

    # ── figure ─────────────────────────────────────────────────────────────
    plt.close('all')
    set_style('science')

    fig = plt.figure(figsize=(7, 6.))
    axd = fig.subplot_mosaic(
        """
        AAAB
        AAAB
        ....
        CDEF
        GHIJ
        GHIJ
        """,
        gridspec_kw={
            'height_ratios': [1, 1, 0.05, 1.8, 1, 1],
            'hspace': 0.55,
            'wspace': 0.5,
        },
    )

    # ── Panel A: flux vs time ───────────────────────────────────────────────
    ax_A = axd['A']
    A_TREF = 3694                        # reference epoch in BTJD
    A_xlim_raw = (3694.25, 3698.5)
    A_xlim = (A_xlim_raw[0] - A_TREF, A_xlim_raw[1] - A_TREF)
    sel = (x_obs >= A_xlim_raw[0]) & (x_obs <= A_xlim_raw[1])
    y_pct = 1e2 * (y_obs - np.nanmedian(y_obs))   # Δ Flux [%]

    # gray raw scatter
    ax_A.scatter(x_obs[sel] - A_TREF, y_pct[sel],
                 s=0.5, c='gray', linewidths=0, rasterized=True, zorder=1)

    # black 20-min time-binned points connected by a thin line
    dt_bin = 20.0 / (24.0 * 60.0)       # 30 min in days
    _t = x_obs[sel]
    _f = y_pct[sel]
    bin_idx = np.floor((_t - _t.min()) / dt_bin).astype(int)
    t_bin, f_bin = [], []
    for bi in np.unique(bin_idx):
        m = bin_idx == bi
        if m.sum() >= 2:
            t_bin.append(_t[m].mean())
            f_bin.append(_f[m].mean())
    t_bin = np.array(t_bin) - A_TREF
    f_bin = np.array(f_bin)
    fin = np.isfinite(t_bin) & np.isfinite(f_bin)
    t_bin, f_bin = t_bin[fin], f_bin[fin]
    ax_A.scatter(t_bin, f_bin, s=3, c='k', linewidths=0, zorder=3)
    # draw connecting line segment-by-segment, skipping gaps > 0.1 days
    gap_mask = np.concatenate(([True], np.diff(t_bin) > 0.1))
    group_ids = np.cumsum(gap_mask)
    for gid in np.unique(group_ids):
        g = group_ids == gid
        ax_A.plot(t_bin[g], f_bin[g], '-', c='k', lw=0.5, zorder=2)

    ax_A.set_xlim(A_xlim)
    ax_A.set_ylabel(r'$\Delta$ Flux [%]', fontsize='x-small')
    ax_A.set_xlabel(f'Time - 2457{A_TREF} [BTJD days]', fontsize='x-small')
    ax_A.tick_params(labelsize='x-small')

    # green bars for spectral coverage (drawn near the bottom of the flux range)
    bar_height = 0.4                              # bar height in %
    y_bar_lo = np.nanpercentile(y_pct[sel], 1) - bar_height
    y_bar_hi = y_bar_lo + bar_height
    for datestr in DATES:
        tmin = spec[datestr]['spectimes'].min() - A_TREF
        tmax = spec[datestr]['spectimes'].max() - A_TREF
        ax_A.fill_between(
            [tmin, tmax], y_bar_lo, y_bar_hi,
            color='green', linewidth=0, zorder=5,
        )

    # ── Panel B: flux vs phase ─────────────────────────────────────────────
    ax_B = axd['B']
    plot_phased_light_curve(
        lc_times, lc_fluxs, T0, PERIOD, None,
        fig=fig, ax=ax_B,
        binsize_phase=BINSIZE_PHASE,
        xlim=[-0.6, 0.6],
        ylim=YLIM,
        showtext=None,
        showtitle=False,
        savethefigure=False,
        c0='darkgray', alpha0=0.15,
        c1='k', alpha1=1,
        phasewrap=True,
        BINMS=1.5,
    )
    ax_B.set_ylim(YLIM)
    ax_B.set_xlim(-0.6, 0.6)
    ax_B.set_ylabel(r'$\Delta$ Flux [%]', fontsize='x-small')
    ax_B.set_xlabel(r'Phase, $\varphi$', fontsize='x-small')
    ax_B.tick_params(labelsize='x-small')

    # ── Panels C/D/E/F + G/H/I/J ──────────────────────────────────────────
    phase_keys = ['C', 'D', 'E', 'F']   # flux vs phase + linecore
    river_keys = ['G', 'H', 'I', 'J']   # specriver

    colorbars = {}

    for pk, rk, datestr, datelabel in zip(
        phase_keys, river_keys, DATES, DATE_LABELS
    ):
        ax_p = axd[pk]
        ax_r = axd[rk]
        sd = spec[datestr]

        # ── CDEF: flux vs phase + linecore ─────────────────────────────────
        xmin_p, xmax_p, ax2 = _plot_phase_linecore(ax_p, lc_times, lc_fluxs, sd, datestr)

        # date title above the phase panel
        ax_p.set_title(datelabel, fontsize='x-small', pad=3)

        # drop x-labels on all CDEF panels
        ax_p.set_xlabel('')

        # drop primary (left) y-label on D, E, F (keep tick labels)
        if pk != 'C':
            ax_p.set_ylabel('')

        # drop secondary (right) Hα y-label on C, D, E (keep tick labels)
        if pk != 'F':
            ax2.set_ylabel('')

        # ── GHIJ: specriver ────────────────────────────────────────────────
        c_mappable, xmin_r, xmax_r = _plot_specriver(ax_r, fig, sd)
        colorbars[rk] = c_mappable

        # sync xlim of phase panel to specriver
        ax_p.set_xlim(xmin_r, xmax_r)
        # also sync twinx xlim
        for child in ax_p.get_shared_x_axes().get_siblings(ax_p):
            child.set_xlim(xmin_r, xmax_r)

        # drop Δv/v_eq y-label on H, I, J (keep tick labels)
        if rk != 'G':
            ax_r.set_ylabel('')

        # inset colorbar only on J
        if rk == 'J':
            cax = ax_r.inset_axes([1.03, 0.01, 0.04, 0.35])
            cb = fig.colorbar(c_mappable, cax=cax, orientation='vertical', extend='both')
            cb.set_label(r"$f_\lambda$ - $f_{\langle t \rangle}$",
                         rotation=90, labelpad=0.5, fontsize='xx-small')
            cb.set_ticks([VMIN_RIVER, VMAX_RIVER])
            cb.set_ticklabels([str(VMIN_RIVER), str(VMAX_RIVER)])
            cb.ax.tick_params(labelsize='xx-small')
            cb.update_ticks()

    # enforce consistent ylim on panel B (must come after all other axes work)
    axd['B'].set_ylim(YLIM)

    # ── save ───────────────────────────────────────────────────────────────
    outpath = join(PLOTDIR, 'mage_mosaic.png')
    savefig(fig, outpath, dpi=300)
    print(f'Saved {outpath}')

    outpath_pdf = outpath.replace('.png', '.pdf')
    fig.savefig(outpath_pdf, bbox_inches='tight', dpi=300)
    print(f'Saved {outpath_pdf}')


if __name__ == '__main__':
    plot_mage_mosaic()
