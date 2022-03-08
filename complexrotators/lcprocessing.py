"""
Contents:
    cr_periodsearch
"""
import numpy as np, pandas as pd
from numpy import array as nparr

import os, multiprocessing, pickle
from complexrotators.paths import RESULTSDIR, DATADIR

from astropy.io import fits
from astrobase import periodbase, checkplot

nworkers = multiprocessing.cpu_count()

def cr_periodsearch(times, fluxs, starid, outdir):
    """
    Given time and flux, run a period-search for objects expected to be complex
    rotators.

    A few plots and pickle files will be written to `outdir` using the `starid`
    string.

    A dictionary of the results is returned, containing:
        'lsp':lsp, 'fine_lsp':fine_lsp, 'times':times, 'fluxs':fluxs,
        'period':fine_lsp['bestperiod'], 't0':np.nanmin(times), 'outdir':outdir
    """

    assert isinstance(starid, str)

    pklpath = os.path.join(outdir, f"{starid}_cr_periodsearch.pkl")
    if os.path.exists(pklpath):
        print(f"Found {pklpath}, loading and continuing.")
        with open(pklpath, 'rb') as f:
            d = pickle.load(f)
        return d

    sep = 1
    if len(times) > 1e4:
        sep = 10
    if len(times) > 1e5:
        sep = 100

    startp, endp = 0.1, 5
    delta_P = 0.2
    stepsize = 1e-5 # for the fine-tuning

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

    print(42*'.')
    print(f"Standard autofreq period: {lsp['bestperiod']:.7f} d")
    print(f"Fine period: {fine_lsp['bestperiod']:.7f} d")
    print(f"Fine - standard: {fine_lsp['bestperiod']-lsp['bestperiod']:.7f} d")
    print(42*'.')

    outfile = os.path.join(
        outdir, f'{starid}_lombscargle_subset_checkplot.png'
    )
    checkplot.checkplot_png(lsp, times, fluxs, fluxs*1e-4,
                            magsarefluxes=True, phasewrap=True,
                            phasesort=True, phasebin=0.002, minbinelems=7,
                            plotxlim=(-0.8,0.8), plotdpi=200,
                            outfile=outfile, verbose=True)

    d = {
        'lsp':lsp, 'fine_lsp':fine_lsp, 'times':times, 'fluxs':fluxs,
        'period':fine_lsp['bestperiod'], 't0':np.nanmin(times), 'outdir':outdir
        }

    with open(pklpath, 'wb') as f:
        pickle.dump(d, f)
        print(f'Made {pklpath}')

    return d
