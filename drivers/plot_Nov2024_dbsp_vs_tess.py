"""
LP 12-502, November 2024: DBSP Hα EW and TESS broadband flux, phase-folded.

Top panel   : TESS Δ Flux [%], colored by orbital cycle number (plasma)
Bottom panel: DBSP Hα EW [Å], colored by orbital cycle number (plasma)

Dates: 20241104 – 20241108.  TESS clipped to DBSP observational window.
Cycle 0 is anchored to the first DBSP exposure on 2024 Nov 04.

env: cpv  (/Users/luke/local/miniconda3/envs/cpv/bin/python)
"""

import os
from astropy.io import fits
from glob import glob
from os.path import join

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy import units as u

from cdips.utils.lcutils import astropy_utc_time_to_bjd_tdb
from aesthetic.plot import set_style, savefig
from complexrotators.paths import RESULTSDIR

# ── LP 12-502 orbital parameters ─────────────────────────────────────────────
T0     = 2450000 + 1791.12 - 0.3 * 18.5611 / 24   # BJD TDB
PERIOD = 18.5611 / 24                               # days

TARGETID    = 'LP_12-502'
RA          = 16.733175 * u.deg
DEC         = 80.45945  * u.deg
UTCDATESTRS = "20241104,20241105,20241106,20241107,20241108".split(",")

DBSP_REDUX = '/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP_REDUX'
TESSDIR    = '/Users/luke/Dropbox/proj/cpv/data/photometry/tess'


def _t_to_phase(t_bjd):
    """Phase in [0, 1)."""
    phi = (t_bjd - T0) / PERIOD
    return phi - np.floor(phi)


def _cycle_number(t_bjd, cycle_ref):
    """Integer cycle number, zero-normalised so cycle_ref maps to 0."""
    return np.floor((t_bjd - T0) / PERIOD).astype(int) - cycle_ref


def _load_ew(utcdatestr):
    """Return (t_bjd, ews_ang, errs_ang) for one night of DBSP Hα."""
    fitsdir   = join(DBSP_REDUX, f'{utcdatestr}/p200_dbsp_red_A/Science')
    ewdir     = join(DBSP_REDUX, f'{utcdatestr}/p200_dbsp_red_A/Balmer_EWs')
    fitspaths = np.sort(glob(join(fitsdir, f'spec1d_red*{TARGETID}*fits')))
    assert len(fitspaths) > 0, f'No FITS for {TARGETID} {utcdatestr}'

    times, keys = [], []
    for f in fitspaths:
        hl   = fits.open(f)
        mjd  = np.float64(hl[0].header['MJD'])
        texp = (hl[0].header['EXPTIME'] * u.second).to(u.day).value
        times.append(mjd + texp)
        keys.append(os.path.basename(f).replace('.fits', ''))
        hl.close()

    csvpaths = [
        glob(join(ewdir, '*HAlpha*' + key + '*_results.csv'))[0]
        for key in keys
    ]
    dfs        = [pd.read_csv(f) for f in csvpaths]
    fitted_ews = np.array([np.abs(df['EW_mA'].iloc[0]) for df in dfs])
    perr       = np.array([df['Fitted_EW_mA_perr'].iloc[0] for df in dfs])
    merr       = np.array([df['Fitted_EW_mA_merr'].iloc[0] for df in dfs])
    errs       = 0.75 * np.array([merr, perr])   # (2, N) milliangstroms

    times = np.array(times) + 2400000.5           # MJD → JD
    t     = Time(times, format='jd', scale='utc')
    t_bjd = astropy_utc_time_to_bjd_tdb(
        t, RA, DEC, observatory='earthcenter', get_barycorr=0
    )
    return t_bjd, 1e-3 * fitted_ews, 1e-3 * errs


def _load_tess(t_bjd_start, t_bjd_end):
    """Return clipped, normalised (t_bjd, dflux_pct, dflux_err_pct)."""
    fitspath = join(
        TESSDIR,
        'tess2024300212641-s0085-0000000402980664-0282-s_lc.fits'
    )
    hdul     = fits.open(fitspath)
    d        = hdul[1].data
    t_bjd    = d['TIME'] + 2457000        # BTJD → BJD TDB
    flux     = d['PDCSAP_FLUX']
    flux_err = d['PDCSAP_FLUX_ERR']
    med      = np.nanmedian(flux)
    flux     /= med
    flux_err /= med

    sel = (t_bjd >= t_bjd_start) & (t_bjd <= t_bjd_end) & np.isfinite(flux)
    t_cl  = t_bjd[sel]
    f_cl  = flux[sel]
    fe_cl = flux_err[sel]

    dflux     = 100.0 * (f_cl - np.nanmedian(f_cl))
    dflux_err = 100.0 * fe_cl
    return t_cl, dflux, dflux_err


def plot_nov2024_dbsp_vs_tess():

    plt.close('all')
    set_style('science')

    # ── load DBSP EW for all nights ───────────────────────────────────────────
    all_t, all_ew, all_err = [], [], []
    for datestr in UTCDATESTRS:
        t, ew, err = _load_ew(datestr)
        all_t.append(t)
        all_ew.append(ew)
        all_err.append(err)

    # anchor cycle 0 to the first DBSP observation on Nov 04
    t_start   = float(np.min(all_t[0]))
    t_end     = float(np.max(all_t[-1]))
    cycle_ref = int(np.floor((t_start - T0) / PERIOD))

    # ── load and clip TESS ────────────────────────────────────────────────────
    t_tess, dflux_tess, derr_tess = _load_tess(t_start, t_end)
    phase_tess = _t_to_phase(t_tess)
    cycle_tess = _cycle_number(t_tess, cycle_ref)

    # simultaneous mask: TESS point falls within any DBSP night's baseline
    sim_mask = np.zeros(len(t_tess), dtype=bool)
    for t_night in all_t:
        sim_mask |= (t_tess >= t_night.min()) & (t_tess <= t_night.max())

    # ── colormap: discrete plasma, one color per cycle ────────────────────────
    vmax     = int(cycle_tess.max())
    n_cycles = vmax + 1
    cmap     = mcolors.ListedColormap(cm.plasma(np.linspace(0, 1, n_cycles)))
    bounds   = np.arange(-0.5, n_cycles + 0.5, 1)
    norm     = mcolors.BoundaryNorm(bounds, cmap.N)

    # ── figure ────────────────────────────────────────────────────────────────
    fig, axs = plt.subplots(
        nrows=2, ncols=1,
        figsize=(3.5, 4.25),
        sharex=True,
    )
    ax_fl, ax_ew = axs

    # top: TESS broadband flux
    # non-simultaneous: small grey semi-transparent background points
    ax_fl.scatter(
        phase_tess[~sim_mask], dflux_tess[~sim_mask],
        c='gray', s=0.8, linewidths=0, zorder=1, alpha=0.25, rasterized=True,
    )
    # simultaneous: colored by cycle number
    ax_fl.scatter(
        phase_tess[sim_mask], dflux_tess[sim_mask],
        c=cycle_tess[sim_mask], cmap=cmap, norm=norm,
        s=1.5, linewidths=0, zorder=2, alpha=0.85, rasterized=True,
    )
    ax_fl.errorbar(
        phase_tess[sim_mask], dflux_tess[sim_mask], yerr=derr_tess[sim_mask],
        fmt='none', ecolor='gray', elinewidth=0.2,
        capsize=0, alpha=0.25, zorder=1,
    )
    ax_fl.set_ylabel(r'$\Delta$ Flux [%]', fontsize='small')
    ax_fl.set_ylim(-2.9, 2.9)

    # bottom: DBSP Hα EW, colored by cycle number (same colormap as TESS)
    for t_bjd, ews, errs in zip(all_t, all_ew, all_err):
        phase = _t_to_phase(t_bjd)
        cycle = _cycle_number(t_bjd, cycle_ref)
        ax_ew.scatter(
            phase, ews,
            c=cycle, cmap=cmap, norm=norm,
            s=10, linewidths=0, zorder=2, alpha=0.9,
        )
        ax_ew.errorbar(
            phase, ews, yerr=errs,
            fmt='none', ecolor='gray', elinewidth=0.5,
            capsize=0, alpha=0.15, zorder=1,
        )
    ax_ew.set_ylabel(r'H$\alpha$ EW [$\mathrm{\AA}$]', fontsize='small')
    ax_ew.set_xlabel(r'Phase ($P$=18.6 hr)', fontsize='small')
    ax_ew.set_xlim(-0.03, 1.03)

    for ax in axs:
        ax.tick_params(axis='both', which='major', labelsize='x-small')

    # leave right margin for colorbar
    fig.tight_layout(rect=[0, 0, 0.84, 1], h_pad=0.4)

    # ── colorbar: right side, half the subplot span, vertically centred ───────
    pos0      = ax_fl.get_position()
    pos1      = ax_ew.get_position()
    span_top  = pos0.y0 + pos0.height
    span_bot  = pos1.y0
    cb_height = (span_top - span_bot) * 0.50
    cb_center = (span_top + span_bot) / 2.0
    cb_y0     = cb_center - cb_height / 2.0
    cb_x0     = pos0.x0 + pos0.width + 0.045
    cb_width  = 0.030

    cbar_ax = fig.add_axes([cb_x0, cb_y0, cb_width, cb_height])
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Cycle', fontsize='x-small')
    cbar.set_ticks(np.arange(0, vmax + 1, 1))
    cbar.ax.tick_params(labelsize='xx-small')

    # ── save ──────────────────────────────────────────────────────────────────
    outdir = join(RESULTSDIR, 'EW_results')
    os.makedirs(outdir, exist_ok=True)
    outpath = join(outdir, 'LP12-502_Nov2024_dbsp_vs_tess.png')
    savefig(fig, outpath, dpi=300)
    print(f'Saved {outpath}')


if __name__ == '__main__':
    plot_nov2024_dbsp_vs_tess()
