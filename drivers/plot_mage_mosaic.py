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
C/D/E/F: flux vs phase + Hα line-core EW twinx (continuum-subtracted EW,
         polycont normalization), for dates 20250118 / ... / 20250121.
G/H/I/J: specriver pcolor removeavg (= sixpanel panel E),
         for dates 20250118 / 20250119 / 20250120 / 20250121.

Variants (see plot_mage_mosaic()):
    avg          : show_simult=False
    simult       : show_simult=True
    simult_tworow: show_simult=True, tworow=True; adds a duplicate G/H/I/J-style
                   row (K/L/M/N) below G/H/I/J, identical except its color
                   scale is vmin/vmax = -0.2/+0.2 instead of -0.2/+0.5.
    simult_tworow_fournightavg:
                   like simult_tworow (show_simult=True, tworow=True,
                   fournightavg=True), but the subtracted "average" template is
                   a SINGLE Gaussian-smoothed 25th-percentile spectrum pooled
                   over ALL four nights (each spectrum interpolated onto a
                   common velocity grid), rather than one template per night.
                   Only the specriver panels (G-J, K-N) change; the C-F
                   linecore twinx is computed before any subtraction and is
                   unaffected.

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
LAMBDA0 = 6562.8                     # angstrom
C_KMS = 299792.458                   # km/s
PERIOD = 0.3439207894624283          # days
T0 = 2460693.243837336 - 2457000     # TJD (BTJD), same as sixpanel driver
YLIM = [-8, 5]                       # flux vs phase ylim [%]
ADOPTED_VEQ = 89.65                  # km/s
DELTA_V = 73.0                       # spectrograph velocity resolution [km/s]
PCTILE = 25                          # for removeavg gaussian smoothed percentile
BINSIZE_PHASE = 0.005
VMIN_RIVER, VMAX_RIVER = -0.2, 0.5  # specriver colormap range
VMIN_RIVER2, VMAX_RIVER2 = -0.2, 0.2  # specriver colormap range, "tworow" duplicate row

# RV sinusoid overlay: v(φ)/v_eq = K_OVER_VEQ * sin(2π(φ − RV_PHI0)).
# The model period IS the phase-fold period PERIOD (0.3439207894624283 d ≈
# 8.254 hr, not the rounded 8.26 hr). Because the x-axis is already phase in
# units of PERIOD, the sinusoid is exactly period-1 in phase, so it tracks
# whatever precise period these plots are folded at — no rescaling needed.
K_OVER_VEQ = 4.0                     # RV semi-amplitude in units of v_eq
RV_PHI0 = -0.12                      # reference phase of the RV sinusoid

DATES = ['20250118', '20250119', '20250120', '20250121']
DATE_LABELS = ['2025/01/18', '2025/01/19', '2025/01/20', '2025/01/21']

LCDIR = join(os.path.expanduser('~'), 'Dropbox/proj/cpv/data/photometry/tess')
LCPATH = join(LCDIR, 'tess2025014115807-s0088-0000000300651846-0285-s_lc.fits')
CACHEDIR = join(LOCALDIR, 'cpv_finding', 'spoc2min_debug')

# 20260820: polycont normalization + EW core metric + duplicate-exposure
# rejection (see mage_continuum_normalization_assessment.txt).  Kept separate
# from results/mage_mosaic/ so the pre-fix plots remain for comparison.
PLOTDIR = join(RESULTSDIR, 'mage_mosaic_20260820')
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
        ew          : continuum-subtracted line-core EW [angstrom] (n_spec,)
        ew_err      : 1-sigma EW uncertainty [angstrom], photon noise plus
                      continuum-fit uncertainty in quadrature (n_spec,)
    """
    (specpaths, spectimes, xvals, yvals, yvalsnonorm, norm_flxs,
     norm_flx_errs) = get_specriver_data(
        TICID, LINESTR, dlambda=DLAMBDA, usespectype='mage', datestr=datestr
    )

    # normalize x to Δv/v_eq
    xvals = [np.asarray(xv) / ADOPTED_VEQ for xv in xvals]

    # getters now matches times to kept files by filename, so lengths agree
    spectimes = np.asarray(spectimes)
    n_spec = len(yvals)
    assert len(spectimes) == n_spec

    # phases centred so that the first observation starts near 0
    xval_phase = (spectimes - T0) / PERIOD
    xval_phase -= np.min(np.ceil(xval_phase))

    # build (n_wav × n_spec) flux array
    n_wav = len(xvals[0])
    flux_arr = np.zeros((n_wav, n_spec))
    for ix, yval in enumerate(yvals):
        flux_arr[:, ix] = yval

    # Continuum-subtracted line-core equivalent width, before avg subtraction:
    #   EW_i = Σ (f_rect − 1) Δλ  over |Δv/v_eq| ≤ 1.
    # (yvals are rectified by the polycont fit, so f_rect ≈ f/f_cont.)  Unlike
    # the previous sum of normalized flux, the ~8-pixel continuum pedestal
    # does not ride along, so a continuum error rescales the line flux instead
    # of adding a pedestal.  Errors: photon noise on the summed counts plus
    # the continuum-fit uncertainty, in quadrature.
    ew, ew_err = [], []
    for yval, ynn, xval, nf, nfe in zip(
        yvals, yvalsnonorm, xvals, norm_flxs, norm_flx_errs
    ):
        dlam = np.nanmedian(np.diff(xval)) * ADOPTED_VEQ * LAMBDA0 / C_KMS
        core = np.abs(xval) <= 1
        # yvalsnonorm skips the ragged-array trim in getters; align defensively
        ynn = np.asarray(ynn)[:len(xval)]
        ew.append(np.sum(yval[core] - 1) * dlam)
        core_counts = np.sum(ynn[core])
        sig_phot = dlam * np.sqrt(core_counts) / nf
        sig_cont = nfe * dlam * np.sum(yval[core])
        ew_err.append(np.sqrt(sig_phot**2 + sig_cont**2))
    ew = np.array(ew)
    ew_err = np.array(ew_err)

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
        'flux_arr': flux_arr,       # raw (n_wav × n_spec), before avg subtraction
        'vel': xvals[0],            # this night's Δv/v_eq grid (n_wav,)
        'flux_arr_sub': flux_arr_sub,
        'norm_flxs': norm_flxs,
        'norm_flx_errs': norm_flx_errs,
        'ew': ew,
        'ew_err': ew_err,
    }


def _apply_fournight_template(spec):
    """Replace each night's per-night 25th-percentile template subtraction with
    a single template pooled over all four nights.

    Because the four nights have slightly different velocity grids (differing
    both in pixel count and in edge velocities), the spectra cannot simply be
    stacked. Instead every spectrum from every night is interpolated onto a
    common velocity grid spanning all nights, the Gaussian-smoothed PCTILE-th
    percentile spectrum is computed once over that pooled set, and it is then
    subtracted from each night after being interpolated back onto that night's
    own native grid.

    Mutates `spec` in place: only 'flux_arr_sub' (used by the specriver panels)
    is overwritten. The linecore/twinx quantities are derived from raw flux and
    are intentionally left untouched.
    """
    # common velocity grid spanning the union of all nights' ranges
    all_vel = [spec[d]['vel'] for d in DATES]
    vmin = min(v.min() for v in all_vel)
    vmax = max(v.max() for v in all_vel)
    n_common = max(len(v) for v in all_vel)
    common_vel = np.linspace(vmin, vmax, n_common)

    # pool every spectrum (all nights) onto the common grid → (n_common, n_tot)
    pooled = []
    for d in DATES:
        v = spec[d]['vel']
        flux_arr = spec[d]['flux_arr']       # raw (n_wav × n_spec)
        for j in range(flux_arr.shape[1]):
            pooled.append(np.interp(common_vel, v, flux_arr[:, j]))
    pooled = np.array(pooled).T

    pctflux = np.nanpercentile(pooled, PCTILE, axis=1)
    common_template = gaussian_filter1d(pctflux, sigma=1)  # sigma=1 for MaGE

    # subtract the pooled template from each night on its own native grid
    for d in DATES:
        sd = spec[d]
        template = np.interp(sd['vel'], common_vel, common_template)
        sd['flux_arr_sub'] = sd['flux_arr'] - template[:, None]

    return spec


def _plot_phase_linecore(ax, lc_times, lc_fluxs, sd, datestr, show_simult=False):
    """Draw flux vs phase + Hα linecore twinx (sixpanel panel D) on ax.

    Replicates the axd['D'] + ax2 linecore logic from
    plot_movie_sixpanel_specriver, style='science' (black on white).

    Returns the specriver phase limits (xmin, xmax) for use in panel G/H/I/J.
    """
    if show_simult:
        # normalize before filtering so the mean is over the full sector
        lc_fluxs = 1e2 * (lc_fluxs - np.nanmean(lc_fluxs))
        t_lo = sd['spectimes'].min()
        t_hi = sd['spectimes'].max()
        mask = (lc_times >= t_lo) & (lc_times <= t_hi)
        lc_times = lc_times[mask]
        lc_fluxs = lc_fluxs[mask]

    # ── linecore twinx (drawn first so TESS LC renders on top) ─────────────
    # continuum-subtracted line-core EW, relative to the night median [%].
    # The 1-σ errors include both photon noise and the continuum-fit
    # uncertainty (no fudge factor).
    ew = sd['ew']
    ew_err = sd['ew_err']

    med_ew = np.nanmedian(ew)
    ew_rel_pct = 100.0 * (ew / med_ew - 1)
    ew_err_pct = 100.0 * ew_err / np.abs(med_ew)

    specphases = sd['specphases']
    xerr = np.nanmedian(np.diff(sd['spectimes']) / PERIOD) / 2

    ax2 = ax.twinx()
    ax2.errorbar(
        specphases, ew_rel_pct,
        yerr=ew_err_pct, xerr=xerr,
        lw=0.8, ls=':', marker='.', c='forestgreen', markersize=2,
        alpha=0.5, zorder=1,
    )
    ax2.set_ylabel(r"$\Delta$ EW$_{\mathrm{H\alpha\ core}}$ [%]",
                   fontsize='x-small', color='forestgreen')
    ax2.tick_params(axis='y', labelcolor='forestgreen', labelsize='x-small')

    # cap outlier-driven scale for 20250119
    #if datestr == '20250119':
    #    ax2.set_ylim(-30, 30)

    # bring TESS LC axis in front of the twinx linecore axis
    ax2.set_zorder(1)
    ax.set_zorder(2)
    ax.patch.set_visible(False)   # keep ax2 background visible

    # ── phased LC ──────────────────────────────────────────────────────────
    if show_simult:
        # Phase in the same coordinate system as specphases
        phase_offset = np.min(np.ceil((sd['spectimes'] - T0) / PERIOD))
        phases = (lc_times - T0) / PERIOD - phase_offset

        # gray raw scatter
        ax.scatter(phases, lc_fluxs, s=0.5, c='gray', linewidths=0,
                   rasterized=True, zorder=1)

        # 20-min time-binned black points (same cadence as panel A)
        dt_bin = 20.0 / (24.0 * 60.0)
        bin_idx = np.floor((lc_times - lc_times.min()) / dt_bin).astype(int)
        t_bin, f_bin = [], []
        for bi in np.unique(bin_idx):
            m = bin_idx == bi
            if m.sum() >= 2:
                t_bin.append(lc_times[m].mean())
                f_bin.append(lc_fluxs[m].mean())
        t_bin = np.array(t_bin)
        f_bin = np.array(f_bin)
        fin = np.isfinite(t_bin) & np.isfinite(f_bin)
        p_bin = (t_bin[fin] - T0) / PERIOD - phase_offset
        ax.scatter(p_bin, f_bin[fin], s=3, c='k', linewidths=0, zorder=3)
        gap_mask = np.concatenate(([True], np.diff(t_bin[fin]) > 0.1))
        group_ids = np.cumsum(gap_mask)
        for gid in np.unique(group_ids):
            g = group_ids == gid
            ax.plot(p_bin[g], f_bin[fin][g], '-', c='k', lw=0.5, zorder=2)
    else:
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


def _plot_specriver(ax, fig, sd, vmin=VMIN_RIVER, vmax=VMAX_RIVER):
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
    delta_v_norm = DELTA_V / ADOPTED_VEQ
    ylim_r = ax_r.get_ylim()
    y_center = ylim_r[0] + 0.14 * (ylim_r[1] - ylim_r[0])
    x_marker = 0.16
    ax_r.errorbar(
        x_marker, y_center, yerr=delta_v_norm / 2,
        transform=trans, fmt='none', ecolor='k',
        elinewidth=0.7, capsize=2, capthick=0.7,
        clip_on=False, zorder=10,
    )
    ax_r.text(
        x_marker + 0.06, y_center, r'$\delta v$',
        transform=trans, rotation=90, ha='left', va='center',
        fontsize='x-small', clip_on=False,
    )


def _overplot_rv_sinusoid(ax):
    """Overplot a thin, semi-transparent velocity sinusoid on a specriver
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


def plot_mage_mosaic(show_simult=False, tworow=False, fournightavg=False):

    # ── load LC ────────────────────────────────────────────────────────────
    (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend,
     cadence_sec, sector, starid) = prepare_local_lc(TICID, LCPATH, CACHEDIR)

    d_ps = cpv_periodsearch(x_obs, y_flat, starid, CACHEDIR, t0='binmin')
    lc_times = d_ps['times']
    lc_fluxs = d_ps['fluxs']

    # ── load spectral data ─────────────────────────────────────────────────
    print('Loading spectral data...')
    spec = {d: _make_specriver_data(d) for d in DATES}

    # optionally replace per-night avg subtraction with a single template
    # pooled over all four nights.
    if fournightavg:
        _apply_fournight_template(spec)

    # ── figure ─────────────────────────────────────────────────────────────
    plt.close('all')
    set_style('science')

    if tworow:
        fig = plt.figure(figsize=(7, 7.7))
        axd = fig.subplot_mosaic(
            """
            AAAB
            AAAB
            ....
            CDEF
            GHIJ
            GHIJ
            KLMN
            KLMN
            """,
            gridspec_kw={
                'height_ratios': [1, 1, 0.05, 1.8, 1, 1, 1, 1],
                'hspace': 0.55,
                'wspace': 0.5,
            },
        )
    else:
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
        xmin_p, xmax_p, ax2 = _plot_phase_linecore(ax_p, lc_times, lc_fluxs, sd, datestr, show_simult=show_simult)

        # date title above the phase panel
        ax_p.set_title(datelabel, fontsize='x-small', pad=3)

        # drop x-labels and tick labels on all CDEF panels
        ax_p.set_xlabel('')
        ax_p.tick_params(axis='x', labelbottom=False)

        # drop primary (left) y-label on D, E, F (keep tick labels)
        if pk != 'C':
            ax_p.set_ylabel('')

        # drop secondary (right) Hα y-label on C, D, E (keep tick labels)
        if pk != 'F':
            ax2.set_ylabel('')

        # ── GHIJ: specriver ────────────────────────────────────────────────
        c_mappable, xmin_r, xmax_r = _plot_specriver(ax_r, fig, sd)
        colorbars[rk] = c_mappable

        # thin semi-transparent RV sinusoid overlay
        _overplot_rv_sinusoid(ax_r)

        # sync xlim of phase panel to specriver
        ax_p.set_xlim(xmin_r, xmax_r)
        # also sync twinx xlim
        for child in ax_p.get_shared_x_axes().get_siblings(ax_p):
            child.set_xlim(xmin_r, xmax_r)

        # drop Δv/v_eq y-label on H, I, J (keep tick labels)
        if rk != 'G':
            ax_r.set_ylabel('')

        # when a KLMN duplicate row follows, GHIJ is no longer the bottom
        # row: drop its x-label/tick labels (same treatment CDEF gets above)
        if tworow:
            ax_r.set_xlabel('')
            ax_r.tick_params(axis='x', labelbottom=False)

        # inset colorbar only on J
        if rk == 'J':
            _add_specriver_colorbar(fig, ax_r, c_mappable)

        # velocity-resolution errorbar marker, lower-left of the bottom-left
        # specriver panel.
        if rk == 'G':
            _add_velocity_marker(ax_r)

    # ── Panels K/L/M/N: duplicate specriver row, vmin/vmax = -0.2/+0.2 ─────
    if tworow:
        river2_keys = ['K', 'L', 'M', 'N']   # duplicate specriver row

        for rk2, datestr in zip(river2_keys, DATES):
            ax_r2 = axd[rk2]
            sd = spec[datestr]

            c_mappable2, xmin_r2, xmax_r2 = _plot_specriver(
                ax_r2, fig, sd, vmin=VMIN_RIVER2, vmax=VMAX_RIVER2
            )

            # thin semi-transparent RV sinusoid overlay
            _overplot_rv_sinusoid(ax_r2)

            # drop Δv/v_eq y-label on L, M, N (keep tick labels)
            if rk2 != 'K':
                ax_r2.set_ylabel('')

            # inset colorbar only on N
            if rk2 == 'N':
                _add_specriver_colorbar(fig, ax_r2, c_mappable2)

            # velocity-resolution errorbar marker, lower-left of the
            # bottom-left specriver panel.
            if rk2 == 'K':
                _add_velocity_marker(ax_r2)

    # enforce consistent ylim on panel B (must come after all other axes work)
    axd['B'].set_ylim(YLIM)

    # Tighten the CDEF-GHIJ gap (and GHIJ-KLMN gap, if present) to the same
    # small absolute size. gridspec's scalar hspace already makes the
    # *default* (pre-shift) row-to-row gaps uniform, but shifting GHIJ up by
    # a fixed amount also drags its bottom edge up, which *widens* the gap
    # below it by that same amount. So a naive "+const" shift applied to
    # both GHIJ and KLMN cancels out and leaves the GHIJ-KLMN gap at its
    # (larger) default rather than tightened. Instead, measure the true
    # default gaps first, then solve for the shifts that land both
    # transitions on the same TIGHT_GAP.
    TIGHT_GAP = 0.008  # target gap size, figure-fraction units

    gap_cg_default = axd['C'].get_position().y0 - axd['G'].get_position().y1
    if tworow:
        gap_gk_default = axd['G'].get_position().y0 - axd['K'].get_position().y1

    d_g = gap_cg_default - TIGHT_GAP
    for key in ['G', 'H', 'I', 'J']:
        pos = axd[key].get_position()
        axd[key].set_position([pos.x0, pos.y0 + d_g, pos.width, pos.height])

    if tworow:
        # K must cancel the widening caused by G's shift (+d_g) on top of
        # its own tightening (-TIGHT_GAP from the default), hence the extra
        # + d_g term relative to a naive same-size shift.
        d_k = gap_gk_default + d_g - TIGHT_GAP
        for key in ['K', 'L', 'M', 'N']:
            pos = axd[key].get_position()
            axd[key].set_position([pos.x0, pos.y0 + d_k, pos.width, pos.height])

    # ── save ───────────────────────────────────────────────────────────────
    if tworow:
        suffix = 'simult_tworow'
        if fournightavg:
            suffix = 'simult_tworow_fournightavg'
    else:
        suffix = 'simult' if show_simult else 'avg'
    outpath = join(PLOTDIR, f'mage_mosaic_{suffix}.png')
    savefig(fig, outpath, dpi=300)
    print(f'Saved {outpath}')

    outpath_pdf = outpath.replace('.png', '.pdf')
    fig.savefig(outpath_pdf, bbox_inches='tight', dpi=300)
    print(f'Saved {outpath_pdf}')


if __name__ == '__main__':
    plot_mage_mosaic(show_simult=False)                    # mage_mosaic_avg.png/pdf
    plot_mage_mosaic(show_simult=True)                      # mage_mosaic_simult.png/pdf
    plot_mage_mosaic(show_simult=True, tworow=True)         # mage_mosaic_simult_tworow.png/pdf
    plot_mage_mosaic(show_simult=True, tworow=True, fournightavg=True)  # mage_mosaic_simult_tworow_fournightavg.png/pdf
