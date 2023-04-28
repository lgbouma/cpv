import os
import numpy as np
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'phase_timegroups_mosaic')
if not os.path.exists(PLOTDIR): os.mkdir(PLOTDIR)

def main():

    ticids = [
        "402980664" # the one and only
    ]

    for ticid in ticids:

        outdir = os.path.join(PLOTDIR, f'TIC_{ticid}')
        if not os.path.exists(outdir): os.mkdir(outdir)

        #for manual_period in np.arange(18.5595, 18.5605, 0.0001):
        for manual_period in [18.5604]:

            cp.plot_phase_timegroups_mosaic(
                outdir,
                ticid=f'TIC_{ticid}',
                lc_cadences='2min',
                binsize_phase=0.005,
                t0=1791.15,
                #t0=1791.19,
                #manual_period=18.559/24,
                manual_period=manual_period/24,
                ylim=[-4.5,3],
                #showtitle=1
                showtitle=0
            )

if __name__ == "__main__":
    main()
