"""
DBSP H-alpha mosaic supplement: 2-column × 2-row overview of EW and
photometric variability for LP 12-502 on 2024 Nov 05 and 2024 Nov 06.

Top row   : Hα EW [Å] vs time
Bottom row: Δ Flux [%] vs time

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

COLUMNS = [
    dict(
        mode='time',
        targetid='LP_12-502',
        utcdatestr='20241105',
        photinst='tess',
        ra=16.733175 * u.deg,
        dec=80.45945 * u.deg,
        label='LP 12-502\n2024 Nov 05',
        xlim_min=-1.2,
    ),
    dict(
        mode='time',
        targetid='LP_12-502',
        utcdatestr='20241106',
        photinst='tess',
        ra=16.733175 * u.deg,
        dec=80.45945 * u.deg,
        label='LP 12-502\n2024 Nov 06',
        xlim_min=-1.2,
    ),
    dict(
        mode='time',
        targetid='LP_12-502',
        utcdatestr='20241108',
        photinst='tess',
        ra=16.733175 * u.deg,
        dec=80.45945 * u.deg,
        label='LP 12-502\n2024 Nov 08',
        xlim_min=-1.2,
    ),
]

DBSP_REDUX = '/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP_REDUX'
TESSDIR    = '/Users/luke/Dropbox/proj/cpv/data/photometry/tess'


def _load_ew_data(targetid, utcdatestr, ra, dec):
    """Load DBSP Hα EW time-series for one night."""
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
    """Load TESS photometry for a single-night time-axis column."""
    assert photinst == 'tess'
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

    t_hours        = 24.0 * (time - 2457000 - x0_btjd)
    delta_flux_pct = 100.0 * (flux - np.nanmedian(flux))
    flux_err_pct   = 100.0 * flux_err

    return t_hours, delta_flux_pct, flux_err_pct, ms


def plot_dbsp_halpha_mosaic_supp():

    plt.close('all')
    set_style('science')

    ncols = len(COLUMNS)
    _F = 0.8
    fig, axs = plt.subplots(
        nrows=2, ncols=ncols,
        figsize=(_F * ncols * 2.0, _F * 4.0),
        sharex='col',
    )

    c = 'k'

    for col_ix, col in enumerate(COLUMNS):
        ax_ew = axs[0, col_ix]
        ax_fl = axs[1, col_ix]

        ax_ew.set_title(col['label'], fontsize='x-small', pad=3)

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

        if 'ew_ymin' in col:
            ax_ew.set_ylim(bottom=col['ew_ymin'])

        ax_ew.text(0.97, 0.05, 'DBSP', transform=ax_ew.transAxes,
                   ha='right', va='bottom', fontsize='xx-small')
        ax_fl.text(0.97, 0.05, 'TESS', transform=ax_fl.transAxes,
                   ha='right', va='bottom', fontsize='xx-small')

        t_pad = 0.04 * (t_hours.max() - t_hours.min())
        xlim_min = col.get('xlim_min', t_hours.min() - t_pad)
        ax_ew.set_xlim(xlim_min, t_hours.max() + t_pad)

    # ── axis labels ──────────────────────────────────────────────────────────
    fs = 'small'
    axs[0, 0].set_ylabel(r'H$\alpha$ EW [$\mathrm{\AA}$]', fontsize=fs)
    axs[1, 0].set_ylabel(r'$\Delta$ Flux [%]', fontsize=fs)

    for col_ix in range(ncols):
        axs[1, col_ix].set_xlabel('Time [hours]', fontsize=fs)

    for ax in axs.flatten():
        ax.tick_params(axis='both', which='major', labelsize='x-small')

    fig.tight_layout()

    outdir = join(RESULTSDIR, 'EW_results')
    os.makedirs(outdir, exist_ok=True)
    outpath = join(outdir, 'DBSP_HAlpha_mosaic_supp.png')
    savefig(fig, outpath, dpi=300)
    print(f'Saved {outpath}')


if __name__ == '__main__':
    plot_dbsp_halpha_mosaic_supp()
