"""
Contents:
    | cpv_periodsearch
"""
import numpy as np, pandas as pd
from numpy import array as nparr

import os, multiprocessing, pickle
from complexrotators.paths import RESULTSDIR, DATADIR

from astropy.io import fits
from astrobase import periodbase, checkplot

from astrobase.lcmath import phase_magseries, phase_bin_magseries


nworkers = multiprocessing.cpu_count()

def cpv_periodsearch(times, fluxs, starid, outdir, t0=None,
                    periodogram_method="pdm"):
    """
    Given time and flux, run a period-search for objects expected to be complex
    rotators.

    A few plots and pickle files will be written to `outdir` using the `starid`
    string.

    Args:

        times (np.ndarray):
            Array of times.

        fluxs (np.ndarray):
            Array of fluxes.

        starid (str):
            Identifier used for cacheing.

        outdir (str):
            Path used for cacheing.

        t0 (None, str, int, or float):
            Epoch at which to phase.  None defaults to t0=1618.  Giving the
            string "binmin" defaults to phase-folding, and taking the
            arg-minimum.  Any int or float will be passed as the manual phase.

        periodogram_method (str):
            "pdm" (phase dispersion minimization) or "ls" (lomb-scargle).

    Returns:

        dict : results

            A dictionary of the results, containing:
                'lsp':periodogram results, 'fine_lsp':fine periodogram results,
                'times':times, 'fluxs':fluxs, 'period':fine_lsp['bestperiod'],
                't0': t0, 'outdir':outdir, 'periodogram_method': ...
            Note that just because the keys are "lsp", the actual method being
            used depends on periodogram_method
    """

    assert isinstance(starid, str)

    pklpath = os.path.join(outdir, f"{starid}_cpv_periodsearch.pkl")
    if os.path.exists(pklpath):
        print(f"Found {pklpath}, loading and continuing.")
        with open(pklpath, 'rb') as f:
            d = pickle.load(f)
        return d

    sep = 1
    # eg., a single TESS sector at 2-minute cadence has 2e4 points.  this cuts
    # it down for period-search purposes to 4e3, which helps the runtime!
    if len(times) > 1e4:
        sep = 5
    # eg., a single TESS sector at 2-minute cadence has 1.2e5 points.
    if len(times) > 1e5:
        sep = 50

    startp, endp = 0.05, 5

    # for the fine-tuning
    delta_P = 0.2
    stepsize = 1e-5

    if periodogram_method == 'ls':

        lsp = periodbase.pgen_lsp(
            times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
            startp=startp, endp=endp, autofreq=True, sigclip=5.0
        )

        fine_lsp = periodbase.pgen_lsp(
            times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
            startp=lsp['bestperiod']-delta_P*lsp['bestperiod'],
            endp=lsp['bestperiod']+delta_P*lsp['bestperiod'],
            autofreq=False, sigclip=5.0, stepsize=stepsize
        )

    elif periodogram_method == 'pdm':

        lsp = periodbase.stellingwerf_pdm(
            times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
            startp=startp, endp=endp, autofreq=True, sigclip=5.0
        )

        fine_lsp = periodbase.stellingwerf_pdm(
            times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
            startp=lsp['bestperiod']-delta_P*lsp['bestperiod'],
            endp=lsp['bestperiod']+delta_P*lsp['bestperiod'],
            autofreq=False, sigclip=5.0, stepsize=stepsize
        )

    print(42*'.')
    print(f"Standard autofreq period: {lsp['bestperiod']:.7f} d")
    print(f"Fine period: {fine_lsp['bestperiod']:.7f} d")
    print(f"Fine - standard: {fine_lsp['bestperiod']-lsp['bestperiod']:.7f} d")
    print(42*'.')

    outfile = os.path.join(
        outdir, f'{starid}_{periodogram_method}_subset_checkplot.png'
    )
    checkplot.checkplot_png(lsp, times, fluxs, fluxs*1e-4,
                            magsarefluxes=True, phasewrap=True,
                            phasesort=True, phasebin=0.002, minbinelems=7,
                            plotxlim=(-0.6,0.6), plotdpi=200,
                            outfile=outfile, verbose=True)

    if t0 is None:
        # default phase
        t0 = 1642.

    elif t0 == 'binmin':
        # bin the phase-fold to 50 points, take the minimum index.

        period = fine_lsp['bestperiod']
        x,y = times, fluxs-np.nanmean(fluxs)
        t0_ini = np.nanmin(x)
        _pd = phase_magseries(x, y, period, t0_ini, wrap=False,
                              sort=False)
        x_fold = _pd['phase']
        y = _pd['mags']
        bs_days = period/50
        orb_bd = phase_bin_magseries(x_fold, y, binsize=bs_days, minbinelems=3)
        min_phase = orb_bd['binnedphases'][np.argmin(orb_bd['binnedmags'])]
        t0 = t0_ini + min_phase*period

    elif isinstance(t0, (int, float)):
        pass

    else:
        raise NotImplementedError

    d = {
        'lsp':lsp, 'fine_lsp':fine_lsp, 'times':times, 'fluxs':fluxs,
        'period':fine_lsp['bestperiod'], 't0': t0, 'outdir':outdir,
        'periodogram_method': periodogram_method
        }

    with open(pklpath, 'wb') as f:
        pickle.dump(d, f)
        print(f'Made {pklpath}')

    return d
