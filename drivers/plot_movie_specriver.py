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
    ]
    linestr = 'Hα' # "Hγ"

    paramdict = {
        # period in hr, t0, flux ylim, sector number, lambdaylim
        '141146667': [0.163762133*24, 3339.9326, [-10, 6], 75,
                      {'Hα':[0.7, 3.8], 'Hγ': [0.7, 7], 'CaH': [0.7, 10]}]
    }
    cb_tickd = {
        'Hα': [1,2,3],
        'Hγ': None,
        'CaH': None
    }

    for ticid in ticids:

        outdir = os.path.join(PLOTDIR, f'TIC_{ticid}')
        if not os.path.exists(outdir): os.mkdir(outdir)

        manual_period = paramdict[ticid][0]
        t0 = paramdict[ticid][1]
        ylim = paramdict[ticid][2]
        sector = paramdict[ticid][3]
        lamylim = paramdict[ticid][4][linestr]

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
                style='science_wob',
                arial_font=1,
                cb_ticks=cb_tickd[linestr]
            )

if __name__ == "__main__":
    main()
