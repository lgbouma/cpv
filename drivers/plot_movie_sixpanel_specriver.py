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

    ##########################################
    # USER CHOICE

    # TIC 141146667: j531, j533
    # ticid = "141146667"

    # LP 12-5602: j546, j547
    # ticid = "402980664"
    # datestr = 'j547'

    ticid = '300651846'
    datestr = '20250102'
    ##########################################

    if ticid == '141146667':
        lcdir = ('/Users/luke/Dropbox/proj/cpv/data/photometry/'+
                 'tess/mastDownload/TESS/'+
                 'tess2024030031500-s0075-0000000141146667-0270-s'
                )
        lcpath = join(
            lcdir, "tess2024030031500-s0075-0000000141146667-0270-s_lc.fits"
        )
    elif ticid == '402980664':
        lcdir = '/Users/luke/Dropbox/proj/cpv/data/photometry/tess'
        if datestr in ['j546', 'j547']:
            lcpath = join(
                lcdir, "tess2024326142117-s0086-0000000402980664-0283-s_lc.fits"
            )
        elif datestr in ['j531', 'j533']:
            lcpath = join(
                lcdir, "tess2023341045131-s0073-0000000402980664-0268-s_lc.fits"
            )
        else:
            raise NotImplementedError
    elif ticid == '300651846':
        lcdir = '/Users/luke/Dropbox/proj/cpv/data/photometry/tess'
        lcpath = join(
            lcdir, 'tess2024353092137-s0087-0000000300651846-0284-s_lc.fits'
        )

    linestr = 'Hα' # "Hγ"
    dlambda = 15

    linestr = 'Hβ' # "Hγ"
    dlambda = 15

    paramdict = {
        # period in hr, t0, flux ylim, sector number, lambdaylim, dlambda
        '141146667': [0.163762133*24,
                      3339.9326,
                      [-10, 6], 75,
                      {'Hα':[0.7, 3.6], 'Hγ': [0.7, 7], 'CaH': [0.7, 10]}, 15],
        '402980664': [18.5611,
                      1791.12,
                      [-4.8,3], 73,
                      {'Hα':[0.7, 2], 'Hγ': [0.7, 7], 'CaH': [0.7, 10]}, 10],
        '300651846': [0.3439207894624283*24,
                      2460693.243837336-2457000,
                      [-8,5], 87,
                      {'Hα':[0.8, 2.5], 'Hβ': [0.7, 6.7], 'Hγ': [0.7, 10]}, dlambda],
    }

    cb_tickd = {
        '141146667': {
            'Hα': [1,2], 'Hγ': None, 'CaH': None
        },
        '402980664': {
            'Hα': [1,2], 'Hγ': None, 'CaH': None
        },
        '300651846': {
            'Hα': [1,1.5], 'Hβ': None, 'Hγ': None, 'CaH': None
        }
    }

    ticids = [ticid]
    for ticid in ticids:

        manual_period = paramdict[ticid][0]
        t0 = paramdict[ticid][1]
        ylim = paramdict[ticid][2]
        sector = paramdict[ticid][3]
        lamylim = paramdict[ticid][4][linestr]
        dlambda = paramdict[ticid][5]

        style = 'science_wob' # "science_wob" or "science"

        outdir = os.path.join(PLOTDIR, f'TIC_{ticid}_{datestr}_{style}')
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
            lognorm=0,
            showsinusoid=0,
            datestr=datestr
        )

if __name__ == "__main__":
    main()
