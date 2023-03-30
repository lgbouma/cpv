import os
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'phase_timegroups')
if not os.path.exists(PLOTDIR): os.mkdir(PLOTDIR)

def main():

    ticids = [
        #"224283342",
        #"146539195"
        "402980664"
    ]

    for ticid in ticids:

        outdir = os.path.join(PLOTDIR, f'TIC_{ticid}')
        if not os.path.exists(outdir): os.mkdir(outdir)

        cp.plot_phase_timegroups(
            outdir,
            ticid=f'TIC_{ticid}',
            lc_cadences='2min',
            binsize_minutes=10,
            t0='binmin',
            ylim=[-100,5],
            yoffset=5.5,
            showtitle=0
        )

if __name__ == "__main__":
    main()
