import os
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'phase')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

def main():

    ticids = [
        #"224283342",
        #"146539195"
        #"402980664" # ID-0
        #"300651846" # northern CVZ
        "405910546" # long period
    ]

    for ticid in ticids:

        outdir = os.path.join(PLOTDIR, f'TIC_{ticid}')
        if not os.path.exists(outdir): os.mkdir(outdir)

        cp.plot_phase(
            outdir,
            ticid=f'TIC_{ticid}',
            lc_cadences='2min',
            ylim=None,
            binsize_minutes=5,
            #manual_period=(8.255/24)
        )

if __name__ == "__main__":
    main()
