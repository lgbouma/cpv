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
        #"300651846"
    ]

    for ticid in ticids:

        outdir = os.path.join(PLOTDIR, f'TIC_{ticid}')
        if not os.path.exists(outdir): os.mkdir(outdir)

        cp.plot_phase_timegroups(
            outdir,
            ticid=f'TIC_{ticid}',
            lc_cadences='2min',
            binsize_phase=0.005,
            #t0='binmin',
            # for 4029
            t0=1791.12,
            ylim=[-100,5],
            yoffset=5.5,
            manual_period=18.5611/24,
            # for 3006 2min
            #ylim=[-190,5],
            #yoffset=6.5,
            # for 3006 QLP
            #lc_cadences='QLP',
            #binsize_minutes=30,
            #ylim=[-450,5],
            #figsize_y=12,
            #manual_period=8.254/24,
            showtitle=0
        )

if __name__ == "__main__":
    main()
