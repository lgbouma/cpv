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
        # "402980664" # the one and only
        "300651846"
    ]

    paramdict = {
        # period in hr, t0, ylim, sector range
        '402980664': [18.5611, 1791.12, [-4.8,3], None, 3],  #1791.15 default?
        '300651846': [8.254, 2170.+12*(8.254/24), [-9.6, 4.9], range(61,70), 5]
    }

    for ticid in ticids:

        outdir = os.path.join(PLOTDIR, f'TIC_{ticid}')
        if not os.path.exists(outdir): os.mkdir(outdir)

        manual_period = paramdict[ticid][0]
        t0 = paramdict[ticid][1]
        ylim = paramdict[ticid][2]
        sector_range = paramdict[ticid][3]
        N_cyclestobin = paramdict[ticid][4]

        for r in [0]:
            cp.plot_movie_phase_timegroups(
                outdir,
                ticid=f'TIC_{ticid}',
                lc_cadences='2min',
                binsize_phase=0.005,
                t0=t0,
                manual_period=manual_period/24,
                ylim=ylim,
                showtitle=0,
                rasterized=r,
                sector_range=sector_range,
                N_cyclestobin=N_cyclestobin
            )

if __name__ == "__main__":
    main()
