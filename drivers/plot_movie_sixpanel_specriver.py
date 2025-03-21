"""
Make movie of CPV flux vs phase, spectrum river, and line cutout.
"""

import os
from os.path import join
import numpy as np
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'movie_sixpanel_specriver')
if not os.path.exists(PLOTDIR): os.mkdir(PLOTDIR)

def main():

    ticids = [
        "141146667"
    ]

    lcdir = ('/Users/luke/Dropbox/proj/cpv/data/photometry/'+
             'tess/mastDownload/TESS/'+
             'tess2024030031500-s0075-0000000141146667-0270-s'
            )
    lcpath = join(
        lcdir, "tess2024030031500-s0075-0000000141146667-0270-s_lc.fits"
    )

    linestr = 'Hα' # "Hγ"

    paramdict = {
        # period in hr, t0, flux ylim, sector number, lambdaylim, dlambda
        '141146667': [0.163762133*24,
                      3339.9326,
                      [-10, 6], 75,
                      {'Hα':[0.7, 3.6], 'Hγ': [0.7, 7], 'CaH': [0.7, 10]}, 15],
    }
    cb_tickd = {
        '141146667': {
            'Hα': [1,2], 'Hγ': None, 'CaH': None
        },
    }

    for ticid in ticids:

        manual_period = paramdict[ticid][0]
        t0 = paramdict[ticid][1]
        ylim = paramdict[ticid][2]
        sector = paramdict[ticid][3]
        lamylim = paramdict[ticid][4][linestr]
        dlambda = paramdict[ticid][5]

        style = 'science' # "science_wob" or "science"

        outdir = os.path.join(PLOTDIR, f'TIC_{ticid}_{style}')
        if not os.path.exists(outdir): os.mkdir(outdir)

        cp.plot_movie_sixpanel_specriver(
            outdir,
            lcpath,
            ticid=ticid,
            linestr=linestr,
            lc_cadences='2min',
            binsize_phase=0.005,
            t0=t0,
            manual_period=manual_period/24,
            ylim=ylim,
            lamylim=lamylim,
            showtitle=0,
            rasterized=0, # rasterize as pdf?  janky cbars
            sector=sector,
            #showhline=0,
            removeavg=1,
            style=style,
            cb_ticks=cb_tickd[ticid][linestr],
            arial_font=1,
            dlambda=dlambda,
            lognorm=0
        )

if __name__ == "__main__":
    main()
