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

    lcpipeline = 'tars'
    lcpipeline = 'spoc2min'

    ticids = [
        #"402980664" # LP 12-502
        #"300651846" # CPV in CZV
        #"141146667" # spectral torus
        #'294328887' # CPV in CVZ
        '268971806' # B star
    ]

    # period in hr, t0, ylim, sector range, N_cyclestobin, binsize_phase, bin_marker_size, smallms
    paramdict = {
        '402980664': [18.5611, 1791.12, [-4.8,3], None, 3, 0.005, 2, 0.3, None],  #1791.15 default?
        # tic3006 spoc2min s61-s70
        #'300651846': [8.254, 2170.+12*(8.254/24), [-9.6, 4.9], range(61,70), 5, 0.005, 2, 0.3, None],
        # tars version
        '300651846': [8.254, 1325.4461206, [-8.6, 4.9], range(1,96), 5, 0.01, 4, 0.5, None],
        '141146667': [3.930, 2420, [-8, 8], None, 3, 0.01, 6, 0.6, 6],
        '294328887': [8.50804, 1325.446, [-8.6, 4.1], range(1,96), 5, 0.01, 4, 0.5, None],
        '268971806': [31.925761, 1492.7986, [-0.2, 0.2], range(1,96), 2, 0.005, 2, 0.3, None]
    }

    for ticid in ticids:

        outdir = os.path.join(PLOTDIR, f'TIC_{ticid}')
        if not os.path.exists(outdir): os.mkdir(outdir)

        manual_period = paramdict[ticid][0]
        t0 = paramdict[ticid][1]
        ylim = paramdict[ticid][2]
        sector_range = paramdict[ticid][3]
        N_cyclestobin = paramdict[ticid][4]
        binsize_phase = paramdict[ticid][5]
        bin_marker_size = paramdict[ticid][6]
        alpha0 = paramdict[ticid][7]
        raw_marker_size = paramdict[ticid][8]

        for r in [0]:
            cp.plot_movie_phase_timegroups(
                outdir,
                ticid=f'TIC_{ticid}',
                lcpipeline=lcpipeline,
                binsize_phase=binsize_phase,
                bin_marker_size=bin_marker_size,
                raw_marker_size=raw_marker_size,
                alpha0=alpha0,
                t0=t0,
                manual_period=manual_period/24,
                ylim=ylim,
                showtitle=0,
                rasterized=r,
                sector_range=sector_range,
                N_cyclestobin=N_cyclestobin,
                style='science_wob',
                arial_font=1
            )

if __name__ == "__main__":
    main()
