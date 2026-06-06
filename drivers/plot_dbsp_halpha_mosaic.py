"""
DBSP H-alpha mosaic: 4-column × 2-row overview of EW and photometric
variability for DG CVn (2024 Apr 28) and LP 12-502 (2024 Nov 04,
2024 Nov 07, 2023 Nov 11+12).

Top row   : Hα EW [Å] vs time (cols 0-2) or phase (col 3)
Bottom row: Δ Flux [%] vs time (cols 0-2) or phase (col 3)

Style: black-on-white ('science').

env: cpv  (/Users/luke/local/miniconda3/envs/cpv/bin/python)
"""

import os
from astropy.io import fits
from glob import glob
from os.path import join

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy import units as u

from cdips.utils.lcutils import astropy_utc_time_to_bjd_tdb
from aesthetic.plot import set_style, savefig
from complexrotators.paths import RESULTSDIR

# ── LP 12-502 phase parameters (from get_ew_timeseries.t_to_phase) ──────────
T0_LP   = 2450000 + 1791.12 - 0.3 * 18.5611 / 24  # BJD TDB
PERIOD_LP = 18.5611 / 24                            # days

# ── column definitions ───────────────────────────────────────────────────────
# mode='time'  → x-axis is hours from start of night (single date)
# mode='phase' → x-axis is orbital phase (utcdatestrs lists multiple dates)
COLUMNS = [
    dict(
        mode='time',
        targetid='DG_CVn',
        utcdatestr='20240428',
        photinst='tierras',
        ra=202.94307059819 * u.deg,
        dec=29.27621104174 * u.deg,
        label='DG CVn\n2024 Apr 28',
    ),
    dict(
        mode='time',
        targetid='LP_12-502',
        utcdatestr='20241104',
        photinst='tess',
        ra=16.733175 * u.deg,
        dec=80.45945 * u.deg,
        label='LP 12-502\n2024 Nov 04',
    ),
    dict(
        mode='time',
        targetid='LP_12-502',
        utcdatestr='20241107',
        photinst='tess',
        ra=16.733175 * u.deg,
        dec=80.45945 * u.deg,
        label='LP 12-502\n2024 Nov 07',
        ew_ymin=5.45,
    ),
    dict(
        mode='phase',
        targetid='LP_12-502',
        utcdatestrs=['20231111', '20231112'],
        photinst='tierras',
        ra=16.733175 * u.deg,
        dec=80.45945 * u.deg,
        label='LP 12-502\n2023 Nov 11+12 (● +▲)',
        # filled circle for night 1, filled triangle-up for night 2
        markers=['o', '^'],
    ),
]

DBSP_REDUX = '/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP_REDUX'
TIERRASDIR = '/Users/luke/Dropbox/proj/cpv/data/photometry/TIERRAS'
TESSDIR    = '/Users/luke/Dropbox/proj/cpv/data/photometry/tess'


def _t_to_phase_lp(t_bjd):
    """Orbital phase for LP 12-502, wrapped to [0, 1)."""
    phi = (t_bjd - T0_LP) / PERIOD_LP
    return phi - np.floor(phi)


def _load_ew_data(targetid, utcdatestr, ra, dec):
    """Load DBSP Hα EW time-series for one night.

    Returns
    -------
    t_hours  : ndarray        hours since first exposure (BJD TDB zero-point)
    ews_ang  : ndarray        EW in Angstroms
    errs_ang : ndarray (2,N)  asymmetric errors in Angstroms
    x0_btjd  : float          reference epoch (BJD - 2457000)
    t_bjd    : ndarray        absolute BJD TDB times
    """
    fitsdir = join(DBSP_REDUX, f'{utcdatestr}/p200_dbsp_red_A/Science')
    ewdir   = join(DBSP_REDUX, f'{utcdatestr}/p200_dbsp_red_A/Balmer_EWs')
    fitspaths = np.sort(glob(join(fitsdir, f'spec1d_red*{targetid}*fits')))
    assert len(fitspaths) > 0, f'No FITS for {targetid} {utcdatestr}'

    times, keys = [], []
    for f in fitspaths:
        hl = fits.open(f)
        mjd  = np.float64(hl[0].header['MJD'])
        texp = (hl[0].header['EXPTIME'] * u.second).to(u.day).value
        times.append(mjd + texp)
        keys.append(os.path.basename(f).replace('.fits', ''))
        hl.close()

    csvpaths = [
        glob(join(ewdir, '*HAlpha*' + key + '*_results.csv'))[0]
        for key in keys
    ]
    dfs = [pd.read_csv(f) for f in csvpaths]
    fitted_ews = np.array([np.abs(df['EW_mA'].iloc[0]) for df in dfs])
    perr = np.array([df['Fitted_EW_mA_perr'].iloc[0] for df in dfs])
    merr = np.array([df['Fitted_EW_mA_merr'].iloc[0] for df in dfs])
    errs = 0.75 * np.array([merr, perr])   # shape (2, N), milliangstroms

    times = np.array(times) + 2400000.5    # MJD → JD
    t = Time(times, format='jd', scale='utc')
    t_bjd = astropy_utc_time_to_bjd_tdb(
        t, ra, dec, observatory='earthcenter', get_barycorr=0
    )

    x0_btjd = np.nanmin(t_bjd - 2457000)
    t_hours  = 24.0 * (t_bjd - 2457000 - x0_btjd)
    ews_ang  = 1e-3 * fitted_ews
    errs_ang = 1e-3 * errs

    return t_hours, ews_ang, errs_ang, x0_btjd, t_bjd


def _load_phot_data_time(utcdatestr, photinst, x0_btjd):
    """Load photometry for a single-night time-axis column.

    Returns (t_hours, delta_flux_pct, flux_err_pct, ms).
    """
    if photinst == 'tess':
        fitspath = join(
            TESSDIR,
            'tess2024300212641-s0085-0000000402980664-0282-s_lc.fits'
        )
        hdul = fits.open(fitspath)
        d = hdul[1].data
        time     = d['TIME'] + 2457000
        flux     = d['PDCSAP_FLUX']
        flux_err = d['PDCSAP_FLUX_ERR']
        med = np.nanmedian(flux)
        flux     /= med
        flux_err /= med
        ms = 1

    elif photinst == 'tierras':
        if utcdatestr == '20240428':
            df = pd.read_csv(
                join(TIERRASDIR, 'TIC368129164_global_lc.csv'), comment='#'
            )
            sel = (df['BJD TDB'] > 2460428.81) & (df['BJD TDB'] < 2460429.2)
            df = df[sel].rename(
                {'Flux': 'Target Relative Flux',
                 'Flux Error': 'Target Relative Flux Error'},
                axis='columns'
            )
        elif utcdatestr == '20231111':
            df = pd.read_csv(join(
                TIERRASDIR,
                '20231111_TIC402980664_circular_fixed_ap_phot_13.csv'
            ))
        else:
            raise NotImplementedError(f'No Tierras handler for {utcdatestr}')
        assert len(df) > 0
        time     = np.array(df['BJD TDB'])
        flux     = np.array(df['Target Relative Flux'])
        flux_err = np.array(df['Target Relative Flux Error'])
        med = np.nanmedian(flux)
        flux     /= med
        flux_err /= med
        ms = 1

    else:
        raise NotImplementedError(photinst)

    t_hours        = 24.0 * (time - 2457000 - x0_btjd)
    delta_flux_pct = 100.0 * (flux - np.nanmedian(flux))
    flux_err_pct   = 100.0 * flux_err

    return t_hours, delta_flux_pct, flux_err_pct, ms


def _load_phot_tierras_phase(utcdatestr):
    """Load Tierras photometry for a phase-axis column.

    Returns (t_bjd, delta_flux_pct, flux_err_pct) where delta_flux_pct is
    100*(flux/median - 1) so each night is independently centred on 0.
    """
    fns = {
        '20231111': '20231111_TIC402980664_circular_fixed_ap_phot_13.csv',
        '20231112': '20231112_TIC402980664_circular_fixed_ap_phot_21.csv',
    }
    df       = pd.read_csv(join(TIERRASDIR, fns[utcdatestr]))
    t_bjd    = np.array(df['BJD TDB'])
    flux     = np.array(df['Target Relative Flux'])
    flux_err = np.array(df['Target Relative Flux Error'])
    med      = np.nanmedian(flux)
    pct      = 100.0 * (flux / med - 1.0)
    pct_err  = 100.0 * flux_err / med
    return t_bjd, pct, pct_err


def plot_dbsp_halpha_mosaic():

    plt.close('all')
    set_style('science')

    ncols = len(COLUMNS)
    _F = 0.8
    fig, axs = plt.subplots(
        nrows=2, ncols=ncols,
        figsize=(_F*ncols * 2.0, _F*4.0),
        sharex='col',
    )

    c = 'k'

    for col_ix, col in enumerate(COLUMNS):
        ax_ew = axs[0, col_ix]
        ax_fl = axs[1, col_ix]

        # ── title ──────────────────────────────────────────────────────────
        ax_ew.set_title(col['label'], fontsize='x-small', pad=3)

        if col['mode'] == 'time':
            # ── single-night, time x-axis ─────────────────────────────────
            t_hours, ews_ang, errs_ang, x0_btjd, _ = _load_ew_data(
                col['targetid'], col['utcdatestr'], col['ra'], col['dec']
            )
            ax_ew.errorbar(
                t_hours, ews_ang, yerr=errs_ang,
                fmt='o', color=c, ecolor=c,
                elinewidth=0.5, capsize=0, capthick=0, markersize=1.5
            )

            t_ph, dflux, dflux_err, ms = _load_phot_data_time(
                col['utcdatestr'], col['photinst'], x0_btjd
            )
            ax_fl.errorbar(
                t_ph, dflux, yerr=dflux_err,
                fmt='o', color=c, ecolor=c, mfc=c,
                elinewidth=0.2, capsize=0, capthick=0,
                markersize=ms, mew=0
            )
            ax_fl.set_ylim(-3.1, 3.1)

            # Optional EW y-axis lower bound override
            if 'ew_ymin' in col:
                ax_ew.set_ylim(bottom=col['ew_ymin'])

            # Instrument / photometer labels
            ax_ew.text(0.97, 0.05, 'DBSP', transform=ax_ew.transAxes,
                       ha='right', va='bottom', fontsize='xx-small')
            phot_label = 'TESS' if col['photinst'] == 'tess' else 'Tierras'
            ax_fl.text(0.97, 0.05, phot_label, transform=ax_fl.transAxes,
                       ha='right', va='bottom', fontsize='xx-small')

            # Constrain x to spectroscopic coverage
            t_pad = 0.04 * (t_hours.max() - t_hours.min())
            ax_ew.set_xlim(t_hours.min() - t_pad, t_hours.max() + t_pad)

        else:
            # ── two-night, phase x-axis ────────────────────────────────────
            markers = col['markers']
            N_STITCH = 10   # points used to compute the night-to-night offset

            # Load photometry for both nights and compute stitch offset so
            # that the two nights match in the phase range 0.64–0.68, where
            # both datasets overlap and the transition is most visible.
            phot = {}
            for datestr in col['utcdatestrs']:
                t_p, pct_p, err_p = _load_phot_tierras_phase(datestr)
                si = np.argsort(t_p)
                phot[datestr] = (t_p[si], pct_p[si], err_p[si])

            date0, date1 = col['utcdatestrs']
            t0_p, pct0, err0_p = phot[date0]
            t1_p, pct1, err1_p = phot[date1]

            ph0 = _t_to_phase_lp(t0_p)
            ph1 = _t_to_phase_lp(t1_p)

            STITCH_LO, STITCH_HI = 0.640, 0.680
            valid0 = np.isfinite(pct0) & (np.abs(pct0) < 30)
            valid1 = np.isfinite(pct1) & (np.abs(pct1) < 30)
            m0 = (ph0 >= STITCH_LO) & (ph0 <= STITCH_HI) & valid0
            m1 = (ph1 >= STITCH_LO) & (ph1 <= STITCH_HI) & valid1
            stitch_offset = np.nanmedian(pct0[m0]) - np.nanmedian(pct1[m1])

            # Apply stitch offset to night 2, then shift both nights up
            FLUX_DY = 0.8
            phot[date0] = (t0_p, pct0 + FLUX_DY, err0_p)
            phot[date1] = (t1_p, pct1 + stitch_offset + FLUX_DY, err1_p)

            for datestr, mk in zip(col['utcdatestrs'], markers):

                # EW vs phase
                _, ews_ang, errs_ang, _, t_bjd = _load_ew_data(
                    col['targetid'], datestr, col['ra'], col['dec']
                )
                phase_ew = _t_to_phase_lp(t_bjd)

                ax_ew.scatter(
                    phase_ew, ews_ang, c=c, s=4, marker=mk,
                    linewidths=0, zorder=2, alpha=0.9
                )
                ax_ew.errorbar(
                    phase_ew, ews_ang, yerr=errs_ang,
                    fmt='none', ecolor=c, elinewidth=0.5,
                    capsize=0, alpha=0.5, zorder=1
                )

                # Broadband flux vs phase (stitched)
                t_bjd_ph, pct, pct_err = phot[datestr]
                phase_fl = _t_to_phase_lp(t_bjd_ph)

                ax_fl.scatter(
                    phase_fl, pct, c=c, s=1.5, marker=mk,
                    linewidths=0, zorder=2, alpha=0.9
                )
                ax_fl.errorbar(
                    phase_fl, pct, yerr=pct_err,
                    fmt='none', ecolor=c, elinewidth=0.2,
                    capsize=0, alpha=0.3, zorder=1
                )

            ax_fl.set_ylim(-4.5, 3.5)
            ax_ew.set_xlim(0.18, 1.02)

            # Instrument / photometer labels
            ax_ew.text(0.97, 0.05, 'DBSP', transform=ax_ew.transAxes,
                       ha='right', va='bottom', fontsize='xx-small')
            ax_fl.text(0.97, 0.05, 'Tierras', transform=ax_fl.transAxes,
                       ha='right', va='bottom', fontsize='xx-small')

    # ── axis labels ──────────────────────────────────────────────────────────
    fs = 'small'
    axs[0, 0].set_ylabel(r'H$\alpha$ EW [$\mathrm{\AA}$]', fontsize=fs)
    axs[1, 0].set_ylabel(r'$\Delta$ Flux [%]', fontsize=fs)

    for col_ix, col in enumerate(COLUMNS):
        if col['mode'] == 'time':
            axs[1, col_ix].set_xlabel('Time [hours]', fontsize=fs)
        else:
            axs[1, col_ix].set_xlabel(r'Phase ($P$=18.6 hr)', fontsize=fs)

    for ax in axs.flatten():
        ax.tick_params(axis='both', which='major', labelsize='x-small')

    fig.tight_layout()

    outdir = join(RESULTSDIR, 'EW_results')
    os.makedirs(outdir, exist_ok=True)
    outpath = join(outdir, 'DBSP_HAlpha_mosaic.png')
    savefig(fig, outpath, dpi=300)
    print(f'Saved {outpath}')


if __name__ == '__main__':
    plot_dbsp_halpha_mosaic()
