"""
Make movie of CPV flux vs phase, spectrum river, and line cutout.
"""

import os
import numpy as np
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'movie_specriver')
if not os.path.exists(PLOTDIR): os.mkdir(PLOTDIR)

def main():

    ticids = [
        "141146667"
        #"402980664"
    ]
    linestr = 'Hα' # "Hγ"

    paramdict = {
        # period in hr, t0, flux ylim, sector number, lambdaylim, dlambda
        '141146667': [0.163762133*24, 3339.9326, [-10, 6], 75,
                      {'Hα':[0.7, 2.1], 'Hγ': [0.7, 7], 'CaH': [0.7, 10]}, 20],
        '402980664': [18.5611, 1791.12, [-4.8,3], 73,
                      {'Hα':[0.7, 2], 'Hγ': [0.7, 7], 'CaH': [0.7, 10]}, 10]
    }
    cb_tickd = {
        '141146667': {
            'Hα': [1,2], 'Hγ': None, 'CaH': None
        },
        '402980664': {
            'Hα': [1,2], 'Hγ': None, 'CaH': None
        }
    }

    # NOTE: depends on star.  linear norm, tuning vmin/vmax in "lambdaylim"
    # above preferred
    lognorm = 0

    for ticid in ticids:

        outdir = os.path.join(PLOTDIR, f'TIC_{ticid}')
        if not os.path.exists(outdir): os.mkdir(outdir)

        manual_period = paramdict[ticid][0]
        t0 = paramdict[ticid][1]
        ylim = paramdict[ticid][2]
        sector = paramdict[ticid][3]
        lamylim = paramdict[ticid][4][linestr]
        dlambda = paramdict[ticid][5]

        for r in [0]:
            cp.plot_movie_specriver(
                outdir,
                ticid=ticid,
                linestr=linestr,
                lc_cadences='2min',
                binsize_phase=0.005,
                t0=t0,
                manual_period=manual_period/24,
                ylim=ylim,
                lamylim=lamylim,
                showtitle=0,
                rasterized=r,
                sector=sector,
                # NOTE: parameters below for F2024 proposals, but look good!
                #  style='science',
                #  figsize=(8,2.5),
                #  savepdf=1,
                #  showhline=0,
                #  cb_ticks=[1,2],
                # NOTE: default parameter set below
                verticallayout=1,
                removeavg=1,
                style='science_wob',
                cb_ticks=cb_tickd[ticid][linestr],
                arial_font=1,
                dlambda=dlambda,
                lognorm=lognorm
            )

if __name__ == "__main__":
    main()
