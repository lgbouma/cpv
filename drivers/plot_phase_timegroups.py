import os
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'phase_timegroups')
if not os.path.exists(PLOTDIR): os.mkdir(PLOTDIR)

def main():

    ticids = [
        #"224283342",
        #"146539195"
        #"402980664"
        "300651846"
    ]

    for ticid in ticids:

        outdir = os.path.join(PLOTDIR, f'TIC_{ticid}')
        if not os.path.exists(outdir): os.mkdir(outdir)

        if ticid == '402980664':
            t0 = 1791.12
            ylim = [-100,5]
            yoffset = 5.5
            manual_period = 18.5611/24
            do4029_resid = 1
            lc_cadences = '2min'
            figsize_y = 7
            cyclewindow = None
            mingap = 3/24

        elif ticid == '300651846':
            # 2min
            yoffset = 6.5
            manual_period = 8.254/24
            t0 = 2170. + 12*manual_period
            figsize_y = 7

            ## for 2min firsthalf
            #lc_cadences = '2min'
            #cyclewindow = [-1,630]
            #ylim = [-105,5]

            # for 2min secondthalf
            lc_cadences = '2min'
            cyclewindow = [2290,2680]
            ylim = [-133,5]
            mingap = 3/24

            ## for all available 2min
            #lc_cadences = '2min'
            #cyclewindow = None
            #ylim = [-240,5]

            # for 3006 QLP
            #lc_cadences='QLP'
            #ylim=[-450,5],
            #figsize_y=12,

            do4029_resid = 0

        else:
            t0 = 'binmin'
            ylim = None
            yoffset = None
            manual_period = None
            lc_cadences = '2min_20sec'
            figsize_y = 7
            cyclewindow = None
            mingap = 3/24

            do4029_resid = 0

        cp.plot_phase_timegroups(
            outdir,
            ticid=f'TIC_{ticid}',
            lc_cadences='2min',
            binsize_phase=0.005,
            t0=t0,
            ylim=ylim,
            yoffset=yoffset,
            manual_period=manual_period,
            showtitle=0,
            do4029_resid=do4029_resid,
            figsize_y=figsize_y,
            cyclewindow=cyclewindow,
            mingap=mingap
        )

if __name__ == "__main__":
    main()
