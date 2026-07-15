"""
MIKE mosaic: TIC 300651846 single-night (2025 January 02) spectroscopic
overview.  Mirrors plot_mage_mosaic.py, but single-column since there is only
one night of MIKE data.

Layout:
    A
    B
    B

A: flux vs phase + Hα linecore twinx (= mage mosaic panels C/D/E/F),
   TESS sector 87.
B: specriver pcolor removeavg, color scale vmin/vmax = -0.2/+0.2
   (= the LOWEST row, K/L/M/N, of mage_mosaic_simult_tworow).

Unlike the mage mosaic, there is no "simult" variant: the MIKE night (BTJD
3678.53-3678.80) falls inside the sector 87 mid-sector downlink gap
(3676.16-3679.80), so no simultaneous TESS photometry exists.

env: cpv  (use /Users/luke/local/miniconda3/envs/cpv/bin/python)
"""

import os
from os.path import join
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.transforms import blended_transform_factory
from scipy.ndimage import gaussian_filter1d

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
DELTA_V = 9.0                        # MIKE velocity resolution [km/s]
PCTILE = 25                          # for removeavg gaussian smoothed percentile
SMOOTH_SIGMA = 5                     # gaussian_filter1d sigma for MIKE
BINSIZE_PHASE = 0.005
VMIN_RIVER, VMAX_RIVER = -0.2, 0.2  # specriver colormap range
VEL_BINFAC = 2                      # bin factor along the velocity axis

# RV sinusoid overlay: v(φ)/v_eq = K_OVER_VEQ * sin(2π(φ − RV_PHI0)).
# The model period IS the phase-fold period PERIOD (0.3439207894624283 d ≈
# 8.254 hr, not the rounded 8.26 hr). Because the x-axis is already phase in
# units of PERIOD, the sinusoid is exactly period-1 in phase, so it tracks
# whatever precise period these plots are folded at — no rescaling needed.
K_OVER_VEQ = 4.0                    # RV semi-amplitude in units of v_eq
RV_PHI0 = -0.05                     # reference phase of the RV sinusoid

DATESTR = '20250102'
DATE_LABEL = '2025/01/02 MIKE'

LCDIR = join(os.path.expanduser('~'), 'Dropbox/proj/cpv/data/photometry/tess')
LCPATH = join(LCDIR, 'tess2024353092137-s0087-0000000300651846-0284-s_lc.fits')
CACHEDIR = join(LOCALDIR, 'cpv_finding', 'spoc2min_debug')

PLOTDIR = join(RESULTSDIR, 'mike_mosaic')
os.makedirs(PLOTDIR, exist_ok=True)


def _make_specriver_data(datestr):
    """Load MIKE spectra for one night and return processed arrays.

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
        TICID, LINESTR, dlambda=DLAMBDA, usespectype='mike', datestr=datestr
    )

    # normalize x to Δv/v_eq
    xvals = [np.asarray(xv) / ADOPTED_VEQ for xv in xvals]

    # spectimes may have one extra leading entry vs yvals
    spectimes = np.asarray(spectimes)
    n_spec = len(yvals)
    if len(spectimes) > n_spec:
        spectimes = spectimes[len(spectimes) - n_spec:]

    # phases centred so that the first observation starts near 0; the +1
    # shift puts this night at 0.2-1.0 rather than -0.8-0.0
    xval_phase = (spectimes - T0) / PERIOD
    xval_phase -= np.min(np.ceil(xval_phase))
    xval_phase += 1.0

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
    smoothmeanflux = gaussian_filter1d(pctflux, sigma=SMOOTH_SIGMA)

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


def _plot_phase_linecore(ax, lc_times, lc_fluxs, sd):
    """Draw flux vs phase + Hα linecore twinx on ax.

    Replicates the mage mosaic C/D/E/F panel logic, style='science'.

    Returns the phase limits (xmin, xmax) and the twinx axis.
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


def _bin_velocity(y_vel, flux_arr, binfac=VEL_BINFAC):
    """Block-average the velocity axis by an integer factor.

    y_vel: (n_wav,) velocity grid; flux_arr: (n_wav, n_spec).  Trailing
    pixels are trimmed if n_wav is not divisible by binfac.
    """
    y_vel = np.asarray(y_vel)
    n_trim = (len(y_vel) // binfac) * binfac
    y_binned = y_vel[:n_trim].reshape(-1, binfac).mean(axis=1)
    f_binned = np.nanmean(
        flux_arr[:n_trim].reshape(-1, binfac, flux_arr.shape[1]), axis=1
    )
    return y_binned, f_binned


def _plot_specriver(ax, fig, sd, vmin=VMIN_RIVER, vmax=VMAX_RIVER):
    """Draw specriver pcolor removeavg on ax.

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

    y_vel, flux_arr_sub = _bin_velocity(y_vel, flux_arr_sub)

    norm = colors.Normalize(vmin=vmin, vmax=vmax)
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


def _add_specriver_colorbar(fig, ax_r, c_mappable):
    """Inset colorbar for a specriver panel, ticked at the mappable's own
    vmin/vmax (so it works for any color-scale range)."""
    vmin, vmax = c_mappable.norm.vmin, c_mappable.norm.vmax
    cax = ax_r.inset_axes([1.03, 0.01, 0.04, 0.35])
    cb = fig.colorbar(c_mappable, cax=cax, orientation='vertical', extend='both')
    cb.set_label(r"$f_\lambda$ - $f_{\langle t \rangle}$",
                 rotation=90, labelpad=0.5, fontsize='xx-small')
    cb.set_ticks([vmin, vmax])
    cb.set_ticklabels([str(vmin), str(vmax)])
    cb.ax.tick_params(labelsize='xx-small')
    cb.update_ticks()


def _add_velocity_marker(ax_r):
    """Velocity-resolution errorbar marker, lower-left of a specriver panel.
    Horizontal position in axes-fraction x; vertical extent is δv/v_eq in the
    panel's own Δv/v_eq data units so it can be read against the specriver.
    """
    trans = blended_transform_factory(ax_r.transAxes, ax_r.transData)
    delta_v_norm = VEL_BINFAC * DELTA_V / ADOPTED_VEQ
    ylim_r = ax_r.get_ylim()
    y_center = ylim_r[0] + 0.14 * (ylim_r[1] - ylim_r[0])
    x_marker = 0.16
    ax_r.errorbar(
        x_marker, y_center, yerr=delta_v_norm / 2,
        transform=trans, fmt='none', ecolor='k',
        elinewidth=0.5, capsize=2, capthick=0.5,
        clip_on=False, zorder=10,
    )
    ax_r.text(
        x_marker + 0.03, y_center, r'$\delta v$',
        transform=trans, rotation=90, ha='left', va='center',
        fontsize='x-small', clip_on=False,
    )


def _overplot_rv_sinusoid(ax):
    """Overplot a thin, semi-transparent velocity sinusoid on the specriver
    panel: v(φ)/v_eq = K_OVER_VEQ * sin(2π(φ − RV_PHI0)), period = 1 phase
    unit. Preserves the panel's existing x/y limits so the pcolor is unaffected.
    """
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    phi = np.linspace(xlim[0], xlim[1], 600)
    v_over_veq = K_OVER_VEQ * np.sin(2 * np.pi * (phi - RV_PHI0))
    ax.plot(phi, v_over_veq, color='k', lw=0.7, ls=':', alpha=0.18, zorder=6)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def plot_mike_mosaic():

    # ── load LC ────────────────────────────────────────────────────────────
    (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend,
     cadence_sec, sector, starid) = prepare_local_lc(TICID, LCPATH, CACHEDIR)

    d_ps = cpv_periodsearch(x_obs, y_flat, starid, CACHEDIR, t0='binmin')
    lc_times = d_ps['times']
    lc_fluxs = d_ps['fluxs']

    # ── load spectral data ─────────────────────────────────────────────────
    print('Loading spectral data...')
    sd = _make_specriver_data(DATESTR)

    # ── figure ─────────────────────────────────────────────────────────────
    plt.close('all')
    set_style('science')

    fig = plt.figure(figsize=(1.75, 4.2))
    axd = fig.subplot_mosaic(
        """
        A
        B
        B
        """,
        gridspec_kw={
            'height_ratios': [1.8, 1, 1],
            'hspace': 0.55,
        },
    )

    ax_p = axd['A']
    ax_r = axd['B']

    # ── Panel A: flux vs phase + linecore ──────────────────────────────────
    xmin_p, xmax_p, ax2 = _plot_phase_linecore(ax_p, lc_times, lc_fluxs, sd)

    # date title above the phase panel
    ax_p.set_title(DATE_LABEL, fontsize='x-small', pad=3)

    # drop x-labels and tick labels on the phase panel (specriver shares x)
    ax_p.set_xlabel('')
    ax_p.tick_params(axis='x', labelbottom=False)

    # ── Panel B: specriver ─────────────────────────────────────────────────
    c_mappable, xmin_r, xmax_r = _plot_specriver(ax_r, fig, sd)

    # thin semi-transparent RV sinusoid overlay
    _overplot_rv_sinusoid(ax_r)

    # sync xlim of phase panel to specriver
    ax_p.set_xlim(xmin_r, xmax_r)
    for child in ax_p.get_shared_x_axes().get_siblings(ax_p):
        child.set_xlim(xmin_r, xmax_r)

    xticks = [0.4, 0.6, 0.8, 1.0]
    ax_r.set_xticks(xticks)
    ax_p.set_xticks(xticks)

    _add_specriver_colorbar(fig, ax_r, c_mappable)
    _add_velocity_marker(ax_r)

    # Tighten the A-B gap to the same small absolute size used in the mage
    # mosaic: shift B up so the gap equals TIGHT_GAP.
    TIGHT_GAP = 0.008  # target gap size, figure-fraction units
    gap_default = ax_p.get_position().y0 - ax_r.get_position().y1
    d_b = gap_default - TIGHT_GAP
    pos = ax_r.get_position()
    ax_r.set_position([pos.x0, pos.y0 + d_b, pos.width, pos.height])

    # ── save ───────────────────────────────────────────────────────────────
    outpath = join(PLOTDIR, 'mike_mosaic.png')
    savefig(fig, outpath, dpi=300)
    print(f'Saved {outpath}')

    outpath_pdf = outpath.replace('.png', '.pdf')
    fig.savefig(outpath_pdf, bbox_inches='tight', dpi=300)
    print(f'Saved {outpath_pdf}')


if __name__ == '__main__':
    plot_mike_mosaic()
