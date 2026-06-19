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

    #lcpipeline = 'tars'
    lcpipeline = 'qlp'
    #lcpipeline = 'spoc2min'

    ticids = [
        #"402980664" # LP 12-502
        #"300651846" # CPV in CZV
        #"368129164" # DG CVn
        "262400835"
        #"141146667" # spectral torus
        #'294328887' # CPV in CVZ
        #'268971806' # B star
        #'55664696'
        #'220476785'
        #'350519637'
    ]

    # period in hr, t0, ylim, sector range, N_cyclestobin,
    # binsize_phase, bin_marker_size, smallms, abbrevname
    paramdict = {
        '402980664': [18.5611, 1791.12, [-4.8,3], None, 3, 0.005, 2,
                      0.3, None, "LP 12-502"],  #1791.15 default?
        # tic3006 spoc2min s61-s70
        #'300651846': [8.254, 2170.+12*(8.254/24), [-9.6, 4.9], range(61,70), 5, 0.005, 2, 0.3, None],
        # tic3006 tars version
        '300651846': [8.254, 1325.4461206, [-8.6, 4.9], range(1,96),
                      5, 0.01, 4, 0.5, None, "TIC 300651846"],
        '368129164': [6.4378, 1930.46, [-1.3, 1.3], range(1,96),
                      3, 0.01, 4, 0.5, None, "DG CVn"],
        '262400835': [7.1573, 1438.2, [-6, 6], range(1,96),
                      5, 0.01, 4, 0.5, None, "TIC 262400835"],
        '141146667': [3.930, 2420, [-8, 8], None, 3, 0.01, 6, 0.6, 6,
                     None],
        '294328887': [8.50804, 1325.446, [-8.6, 4.1], range(1,96), 5,
                      0.01, 4, 0.5, None, None],
        '268971806': [31.925761, 1492.7986, [-0.2, 0.2], range(1,96),
                      2, 0.005, 2, 0.3, None, None],
        '55664696': [19.73260, 1324.9951, [-8,6], range(1,96), 3,
                     0.01, 4, 0.5, None, None],
        '220476785': [6.826901315, 1326, [-1.5,1.5], range(1,96), 7,
                      0.01, 4, 0.5, None, None],
        '350519637': [17.8161986, 1326, [-7,7], range(1,96), 5, 0.01,
                      4, 0.5, None, None],
    }

    # ticid -> (lo_hr, hi_hr) period range of a secondary signal to fit and
    # subtract before phase-folding; ticids absent here get no subtraction.
    secondary_period_ranges = {
        '368129164': (2.58, 2.62),  # DG CVn: second period superposed on dips
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
        title_head = (
            f"{paramdict[ticid][9]}, "+r"$P$"+f"={manual_period:.3f} hr, "+
            r"$t_{0}$"+f"={t0:.2f} BTJD"
        )

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
                arial_font=1,
                title_head=title_head,
                secondary_period_range_hr=secondary_period_ranges.get(ticid, None)
            )

if __name__ == "__main__":
    main()
