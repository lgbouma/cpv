"""
Contents (mostly deprecated because each tries to do too much):
    get_complexrot_data
    get_complexrot_twentysec_data
"""
import numpy as np, pandas as pd
from numpy import array as nparr

import os, multiprocessing, pickle
from complexrotators.paths import RESULTSDIR, DATADIR

from astropy.io import fits
from astrobase import periodbase, checkplot

nworkers = multiprocessing.cpu_count()

from cdips_followup.quicklooktools import (
    get_tess_data, explore_flux_lightcurves, make_periodogram
)

def get_complexrot_data(ticid, kicid=None, hardcsv=None):
    """
    ticid: str or None
    kicid: str or None
    """

    if ticid is not None:
        outdir = os.path.join(RESULTSDIR, 'river', f'tic_{ticid}')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        pklpath = os.path.join(outdir, f'tic_{ticid}_lsinfo.pkl')
    else:
        assert isinstance(kicid, str)
        outdir = os.path.join(RESULTSDIR, 'river', f'kic_{kicid}')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        pklpath = os.path.join(outdir, f'kic_{kicid}_lsinfo.pkl')

    if not os.path.exists(pklpath):

        if hardcsv is not None:
            df = pd.read_csv(hardcsv)
            times, fluxs = nparr(df.time), nparr(df.flux)

        elif ticid is not None:
            data = get_tess_data(ticid, outdir=outdir, spoc=1)
            times, fluxs = explore_flux_lightcurves(
                data, ticid, outdir=outdir, get_lc=1, require_quality_zero=0,
                pipeline='spoc'
            )
        else:
            # # Saul Rappaport's space-separated and whitened format.
            # datapath = os.path.join(DATADIR, 'photometry', 'kepler',
            #                         f'{kicid}_minus_3p3.dat')
            # df = pd.read_csv(datapath, delim_whitespace=True,
            #                  names=['int', 'time', 'flux', 'raw_flux', 'smoothfn'])
            # times, fluxs = np.array(df['time']), np.array(df['flux'])
            ## LGB's whitened format
            # datapath = os.path.join(DATADIR, 'photometry', 'kepler',
            #                         f'kic{kicid}_whitened.csv')
            # df = pd.read_csv(datapath)
            # times, fluxs = np.array(df['time']), np.array(df['flux_r1'])

            ## # NOTE: oneoff
            #datapath = os.path.join(DATADIR, 'photometry', 'kepler',
            #                        f'kepler1627.csv')
            #df = pd.read_csv(datapath)
            #times, fluxs = np.array(df['time']), np.array(df['flux'])
            raise NotImplementedError(
                'You need to homogenize Kepler reading if you plan to keep doing it'
            )


        sep = 1
        if len(times) > 1e4:
            sep = 10
        if len(times) > 1e5:
            sep = 100

        startp, endp = 0.1, 5
        delta_P = 0.2
        stepsize = 1e-5 # for the fine-tuning

        # startp, endp = 0.4037144, 0.4037146
        # delta_P = 1e-5
        # stepsize = 1e-5 # for the fine-tuning

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

        if isinstance(ticid, str):
            outfile = os.path.join(
                outdir, f'tic_{ticid}_lombscargle_subset_checkplot.png'
            )
        else:
            outfile = os.path.join(
                outdir, f'kic_{kicid}_lombscargle_subset_checkplot.png'
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

    with open(pklpath, 'rb') as f:
        d = pickle.load(f)

    return d


def get_complexrot_twentysec_data(ticid='262400835', kicid=None):
    """
    ticid: str or None
    """

    if ticid != '262400835':
        raise NotImplementedError('data getter not yet automated')

    if ticid is not None:
        outdir = os.path.join(RESULTSDIR, 'river', f'tic_{ticid}')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        pklpath = os.path.join(outdir, f'tic_{ticid}_lsinfo.pkl')
    else:
        raise NotImplementedError

    if not os.path.exists(pklpath):

        fitspath = os.path.join(
            DATADIR,
            'photometry/tess/20sec/MAST_2021-06-03T1207/TESS/tess2020324010417-s0032-0000000262400835-0200-a_fast/tess2020324010417-s0032-0000000262400835-0200-a_fast-lc.fits'
        )

        hl = fits.open(fitspath)

        data = [hl[1].data]

        times, fluxs = explore_flux_lightcurves(
            data, ticid, outdir=outdir, get_lc=1, require_quality_zero=0,
            pipeline='spoc', detrend='median'
        )

        sep = 1
        if len(times) > 1e4:
            sep = 10
        if len(times) > 1e5:
            sep = 100

        startp, endp = 0.1, 5
        delta_P = 0.2 # for fine-tuning
        stepsize = 2e-6 # for the fine-tuning

        # startp, endp = 0.4037144, 0.4037146
        # delta_P = 1e-5
        # stepsize = 1e-5 # for the fine-tuning

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
            outdir, f'tic_{ticid}_lombscargle_subset_checkplot.png'
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

    with open(pklpath, 'rb') as f:
        d = pickle.load(f)

    return d
