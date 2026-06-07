"""
HIRES mosaic: LP 12-502 three-night spectroscopic overview.

Layout (3 columns × 2 rows):

    ABC
    DEF

A/B/C : Δ Flux [TESS, primary] + Δ f_{Hα core} [secondary],
         simultaneous with spectroscopy, phased for j531 / j546 / j547
D/E/F : Hα specriver pcolor (Hα − Avg.),
         for j531 / j546 / j547

env: cpv  (/Users/luke/local/miniconda3/envs/cpv/bin/python)
"""

import os
from os.path import join
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.ndimage import gaussian_filter1d

from aesthetic.plot import set_style, savefig

from astropy.io import fits as astropy_fits

from complexrotators.paths import RESULTSDIR, LOCALDIR
from complexrotators.getters import get_specriver_data

# ── constants ──────────────────────────────────────────────────────────────
TICID        = '402980664'
LINESTR      = 'Hα'
DLAMBDA      = 10
PERIOD       = 18.5611 / 24.0    # days
T0           = 1791.12            # BTJD
ADOPTED_VEQ  = 24.2               # km/s  (2π R★ / P_rot)
YLIM         = [-2.1, 2.1]       # Δ Flux [%]
PCTILE       = 25
SIGMA_SMOOTH = 5                  # gaussian_filter1d sigma for removeavg (HIRES)
VMIN_RIVER, VMAX_RIVER = -0.2, 0.5
FUDGE        = 2                  # error-bar inflation for Hα linecore

DATES       = ['j546',        'j547',        'j531'       ]
DATE_LABELS = ['2024 Dec 03', '2024 Dec 10', '2023 Nov 23']

LCDIR = join(os.path.expanduser('~'), 'Dropbox/proj/cpv/data/photometry/tess')
LCPATHS = {
    'j531': join(LCDIR, 'tess2023341045131-s0073-0000000402980664-0268-s_lc.fits'),
    'j546': join(LCDIR, 'tess2024326142117-s0086-0000000402980664-0283-s_lc.fits'),
    'j547': join(LCDIR, 'tess2024326142117-s0086-0000000402980664-0283-s_lc.fits'),
}
FLUXKEYS = {
    'j531': 'PDCSAP_FLUX',
    'j546': 'SAP_FLUX',
    'j547': 'SAP_FLUX',
}
# j547 window has all points flagged with scattered-light (4096) and stray-light
# (32768); accept those bits so the photometry appears in the overlay.
QUALITY_BITMASKS = {
    'j531': 0,
    'j546': 0,
    'j547': 4096 | 32768,
}

CACHEDIR = join(LOCALDIR, 'cpv_finding', 'spoc2min_debug')
os.makedirs(CACHEDIR, exist_ok=True)

PLOTDIR = join(RESULTSDIR, 'hires_lp12502_mosaic')
os.makedirs(PLOTDIR, exist_ok=True)


def _load_tess_lc(lcpath, fluxkey, quality_bitmask=0):
    """Load TESS LC directly from FITS; return (x_obs [BTJD], y_obs_pct [%]).

    quality_bitmask: flag bits to tolerate — points where
    (qual & ~quality_bitmask) == 0 are kept.
    """
    hdul = astropy_fits.open(lcpath)
    data = hdul[1].data
    time = data['TIME']
    flux = data[fluxkey]
    qual = data['QUALITY']
    hdul.close()

    sel = (qual & ~quality_bitmask) == 0
    x_obs = time[sel]
    y_obs = flux[sel].astype(float)
    y_obs /= np.nanmedian(y_obs)
    return x_obs, 1e2 * (y_obs - 1.0)


def _make_specriver_data(datestr):
    """Load HIRES spectra for one night; return processed arrays."""
    specpaths, spectimes, xvals, yvals, yvalsnonorm, norm_flxs = get_specriver_data(
        TICID, LINESTR, dlambda=DLAMBDA, usespectype='reduced', datestr=datestr
    )

    xvals = [np.asarray(xv) / ADOPTED_VEQ for xv in xvals]

    spectimes = np.asarray(spectimes)
    n_spec = len(yvals)
    if len(spectimes) > n_spec:
        spectimes = spectimes[len(spectimes) - n_spec:]

    # phases: centre first observation near 0
    xval_phase = (spectimes - T0) / PERIOD
    xval_phase -= np.min(np.ceil(xval_phase))

    n_wav = len(xvals[0])
    flux_arr = np.zeros((n_wav, n_spec))
    for ix, yval in enumerate(yvals):
        flux_arr[:, ix] = yval

    # linecore: sum flux within |Δv/v_eq| ≤ 1
    linecore_sums = np.array([
        np.sum(yval[np.abs(xval) <= 1])
        for yval, xval in zip(yvals, xvals)
    ])

    # removeavg: subtract gaussian-smoothed 25th-percentile spectrum
    pctflux = np.nanpercentile(flux_arr, PCTILE, axis=1)
    smoothmeanflux = gaussian_filter1d(pctflux, sigma=SIGMA_SMOOTH)

    orig_yvals = deepcopy(yvals)
    for yval in yvals:
        yval -= smoothmeanflux
    flux_arr_sub = flux_arr - smoothmeanflux[:, None]

    return {
        'spectimes'    : spectimes,
        'specphases'   : xval_phase,
        'xvals'        : xvals,
        'orig_yvals'   : orig_yvals,
        'flux_arr_sub' : flux_arr_sub,
        'norm_flxs'    : norm_flxs,
        'linecore_sums': linecore_sums,
    }


def _plot_phase_linecore(ax, x_obs, y_obs_pct, sd):
    """Draw Δ Flux [TESS] (primary) + Hα linecore twinx on ax.

    Returns (xmin, xmax, ax2).
    """
    # ── linecore secondary axis (drawn first so TESS renders on top) ───────
    linecore_sums = sd['linecore_sums']
    norm_flxs     = sd['norm_flxs']
    med_normflx   = np.nanmedian(norm_flxs)

    linecore_rel     = linecore_sums / np.nanmedian(linecore_sums)
    linecore_err     = (np.sqrt(linecore_sums * med_normflx)
                        / (linecore_sums * med_normflx))
    linecore_rel_pct = 100.0 * (linecore_rel - np.nanmedian(linecore_rel))
    linecorr_err_pct = 100.0 * linecore_err

    specphases = sd['specphases']
    xerr = np.nanmedian(np.diff(sd['spectimes']) / PERIOD) / 2

    ax2 = ax.twinx()
    ax2.errorbar(
        specphases, linecore_rel_pct,
        yerr=FUDGE * linecorr_err_pct, xerr=xerr,
        lw=0.8, ls=':', marker='.', c='forestgreen', markersize=2,
        alpha=0.5, zorder=1,
    )
    ax2.set_ylabel(
        r"$\Delta$ $f$$_{\mathrm{H\alpha\ core}}$ [%]",
        fontsize='x-small', color='forestgreen',
    )
    ax2.tick_params(axis='y', labelcolor='forestgreen', labelsize='x-small')

    ax2.set_zorder(1)
    ax.set_zorder(2)
    ax.patch.set_visible(False)

    # ── simultaneous TESS LC ───────────────────────────────────────────────
    t_lo = sd['spectimes'].min()
    t_hi = sd['spectimes'].max()
    mask = (x_obs >= t_lo) & (x_obs <= t_hi)
    lc_t = x_obs[mask]
    lc_f = y_obs_pct[mask]

    phase_offset = np.min(np.ceil((sd['spectimes'] - T0) / PERIOD))
    phases = (lc_t - T0) / PERIOD - phase_offset

    if len(phases) > 0:
        ax.scatter(phases, lc_f, s=0.5, c='gray', linewidths=0,
                   rasterized=True, zorder=1)

        dt_bin  = 20.0 / (24.0 * 60.0)
        bin_idx = np.floor((lc_t - lc_t.min()) / dt_bin).astype(int)
        t_bin, f_bin = [], []
        for bi in np.unique(bin_idx):
            m = bin_idx == bi
            if m.sum() >= 2:
                t_bin.append(lc_t[m].mean())
                f_bin.append(lc_f[m].mean())

        if len(t_bin) > 0:
            t_bin = np.array(t_bin)
            f_bin = np.array(f_bin)
            fin   = np.isfinite(t_bin) & np.isfinite(f_bin)
            p_bin = (t_bin[fin] - T0) / PERIOD - phase_offset
            ax.scatter(p_bin, f_bin[fin], s=3, c='k', linewidths=0, zorder=3)
            gap_mask  = np.concatenate(([True], np.diff(t_bin[fin]) > 0.1))
            group_ids = np.cumsum(gap_mask)
            for gid in np.unique(group_ids):
                g = group_ids == gid
                ax.plot(p_bin[g], f_bin[fin][g], '-', c='k', lw=0.5, zorder=2)

    ax.set_ylim(YLIM)
    xmin = specphases.min() - 0.05
    xmax = specphases.max() + 0.05
    ax.set_xlim(xmin, xmax)
    ax.set_ylabel(r"$\Delta$ Flux [%]", fontsize='x-small')
    ax.set_xlabel(r"Phase, $\varphi$",  fontsize='x-small')
    ax.tick_params(labelsize='x-small')

    return xmin, xmax, ax2


def _plot_specriver(ax, fig, sd):
    """Draw Hα specriver pcolor (Hα − Avg.) on ax.

    Returns (pcolor mappable, xmin, xmax).
    """
    y_vel    = sd['xvals'][0]
    x_phase  = sd['specphases']
    flux_sub = sd['flux_arr_sub']

    norm = colors.Normalize(vmin=VMIN_RIVER, vmax=VMAX_RIVER)
    c = ax.pcolor(
        x_phase, y_vel, flux_sub,
        cmap='Spectral_r', norm=norm, shading='auto', rasterized=True,
    )

    ax.set_ylim(-13, 13)
    ax.set_ylabel(r"$\Delta v/v_\mathrm{eq}$", fontsize='x-small')
    ax.set_xlabel(r"Phase, $\varphi$",          fontsize='x-small')
    ax.tick_params(labelsize='x-small')
    ax.text(
        0.04, 0.95, f'{LINESTR} − Avg.',
        ha='left', va='top', transform=ax.transAxes, fontsize='xx-small',
        bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', pad=0.4),
    )

    xmin, xmax = ax.get_xlim()
    return c, xmin, xmax


def plot_hires_mosaic():
    """Produce the HIRES mosaic for LP 12-502 (avg-profile + 3 × 2 grid)."""

    print('Loading spectral data ...')
    spec = {d: _make_specriver_data(d) for d in DATES}

    print('Loading TESS light curves ...')
    _lc_cache = {}
    lc_data   = {}
    for d in DATES:
        lp  = LCPATHS[d]
        fk  = FLUXKEYS[d]
        qm  = QUALITY_BITMASKS[d]
        key = (lp, fk, qm)
        if key not in _lc_cache:
            _lc_cache[key] = _load_tess_lc(lp, fk, quality_bitmask=qm)
        lc_data[d] = _lc_cache[key]

    # median Hα profile across all dates/spectra = f_⟨t⟩
    all_orig    = [yv for d in DATES for yv in spec[d]['orig_yvals']]
    avg_profile = np.nanmedian(np.array(all_orig), axis=0)
    vel_axis    = spec[DATES[0]]['xvals'][0]

    plt.close('all')
    set_style('science')

    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(5.6, 3.6))
    gs  = GridSpec(
        2, 4, figure=fig,
        width_ratios=[0.23, 1, 1, 1],
        height_ratios=[1, 1.5],
        hspace=0.2, wspace=0.45,
        left=-0.05,
    )

    # axes: profile panel (bottom-left only) + top/bottom triplets
    ax_prof  = fig.add_subplot(gs[1, 0])
    axes_top = [fig.add_subplot(gs[0, c + 1]) for c in range(3)]
    axes_bot = [fig.add_subplot(gs[1, c + 1], sharey=ax_prof) for c in range(3)]

    # ── average Hα profile panel ─────────────────────────────────────────
    ax_prof.plot(avg_profile, vel_axis, 'k-', lw=0.7)
    ax_prof.set_ylim(-13, 13)
    ax_prof.tick_params(labelsize='x-small')
    ax_prof.set_xticks([1, 5])
    ax_prof.invert_xaxis()
    ax_prof.tick_params(axis='y', which='both', labelleft=False, labelright=False)
    ax_prof.set_xlabel(r'$f_{\langle t \rangle}$', fontsize='x-small')

    # ── main loop ────────────────────────────────────────────────────────
    for col, (datestr, datelabel) in enumerate(zip(DATES, DATE_LABELS)):
        ax_p = axes_top[col]
        ax_r = axes_bot[col]
        sd   = spec[datestr]
        x_obs, y_obs_pct = lc_data[datestr]

        # top row: phased TESS + Hα linecore
        xmin_p, xmax_p, ax2 = _plot_phase_linecore(ax_p, x_obs, y_obs_pct, sd)
        ax_p.set_title(datelabel, fontsize='x-small', pad=3)
        ax_p.set_xlabel('')
        ax_p.tick_params(axis='x', labelbottom=False)
        if col != 0:
            ax_p.set_ylabel('')
        if col != 2:
            ax2.set_ylabel('')

        # bottom row: specriver
        c_mappable, xmin_r, xmax_r = _plot_specriver(ax_r, fig, sd)

        # sync xlim
        ax_p.set_xlim(xmin_r, xmax_r)
        ax2.set_xlim(xmin_r, xmax_r)

        # ylabel/ticks on first specriver only
        if col == 0:
            ax_r.set_ylabel(r"$\Delta v/v_\mathrm{eq}$", fontsize='x-small',
                            labelpad=-1)
        else:
            ax_r.set_ylabel('')
            ax_r.tick_params(axis='y', labelleft=False)

        # inset colorbar on rightmost column only
        if col == 2:
            cax = ax_r.inset_axes([1.03, 0.01, 0.04, 0.35])
            cb  = fig.colorbar(c_mappable, cax=cax, orientation='vertical',
                               extend='both')
            cb.set_label(r"$f_\lambda - f_{\langle t \rangle}$",
                         rotation=90, labelpad=0.5, fontsize='xx-small')
            cb.set_ticks([VMIN_RIVER, VMAX_RIVER])
            cb.set_ticklabels([str(VMIN_RIVER), str(VMAX_RIVER)])
            cb.ax.tick_params(labelsize='xx-small')
            cb.update_ticks()

    # tighten gap between top and bottom rows (shift entire bottom row up)
    for ax in [ax_prof] + axes_bot:
        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0 + 0.04, pos.width, pos.height])

    outpath = join(PLOTDIR, 'hires_lp12502_mosaic.png')
    savefig(fig, outpath, dpi=300)
    print(f'Saved {outpath}')

    outpath_pdf = outpath.replace('.png', '.pdf')
    fig.savefig(outpath_pdf, bbox_inches='tight', dpi=300)
    print(f'Saved {outpath_pdf}')


if __name__ == '__main__':
    plot_hires_mosaic()
