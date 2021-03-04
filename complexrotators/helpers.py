"""
Contents:
    get_complexrot_data
"""
import numpy as np, pandas as pd

import os, multiprocessing, pickle
from complexrotators.paths import RESULTSDIR

from astrobase import periodbase, checkplot

nworkers = multiprocessing.cpu_count()

from cdips_followup.quicklooktools import (
    get_tess_data, explore_flux_lightcurves, make_periodogram
)

def get_complexrot_data(ticid):

    outdir = os.path.join(RESULTSDIR, 'river', f'tic_{ticid}')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    pklpath = os.path.join(outdir, f'tic_{ticid}_lsinfo.pkl')

    if not os.path.exists(pklpath):

        data = get_tess_data(ticid, outdir=outdir, spoc=1)

        times, fluxs = explore_flux_lightcurves(
            data, ticid, outdir=outdir, get_lc=1
        )

        sep = 1
        if len(times) > 1e4:
            sep = 10
        if len(times) > 1e5:
            sep = 100

        lsp = periodbase.pgen_lsp(
            times[::sep], fluxs[::sep], fluxs[::sep]*1e-4, magsarefluxes=True,
            startp=0.1, endp=5, autofreq=True, sigclip=5.0
        )

        outfile = os.path.join(
            outdir, f'tic_{ticid}_lombscargle_subset_checkplot.png'
        )

        checkplot.checkplot_png(lsp, times, fluxs, fluxs*1e-4,
                                magsarefluxes=True, phasewrap=True,
                                phasesort=True, phasebin=0.002, minbinelems=7,
                                plotxlim=(-0.8,0.8), plotdpi=200,
                                outfile=outfile, verbose=True)

        d = {
            'lsp':lsp, 'times':times, 'fluxs':fluxs,
            'period':lsp['bestperiod'], 't0':np.nanmin(times), 'outdir':outdir
            }
        with open(pklpath, 'wb') as f:
            pickle.dump(d, f)
            print(f'Made {pklpath}')

    with open(pklpath, 'rb') as f:
        d = pickle.load(f)

    return d
