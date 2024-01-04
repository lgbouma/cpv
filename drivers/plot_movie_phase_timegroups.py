"""
Make figures for a movie of CPV evolution
"""

import os
import numpy as np
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'movie_phase_timegroups')
if not os.path.exists(PLOTDIR): os.mkdir(PLOTDIR)

def main():

    ticids = [
        "402980664" # the one and only
    ]

    for ticid in ticids:

        outdir = os.path.join(PLOTDIR, f'TIC_{ticid}')
        if not os.path.exists(outdir): os.mkdir(outdir)

        for manual_period in [18.5611]:
            for r in [0]:
                cp.plot_movie_phase_timegroups(
                    outdir,
                    ticid=f'TIC_{ticid}',
                    lc_cadences='2min',
                    binsize_phase=0.005,
                    #t0=1791.15, # NOTE default
                    t0=1791.12,
                    #t0=1791.19,
                    manual_period=manual_period/24,
                    ylim=[-4.8,3],
                    showtitle=0,
                    rasterized=r
                )

if __name__ == "__main__":
    main()
