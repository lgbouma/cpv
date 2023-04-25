"""
reproduce a figure from max gunther's 2022 paper, but with 30min vs 10min vs
120sec comparison

his figure was with TIC 141146667, P=3.9hr, an extreme case.  this is probably
fine.
"""

import os
import complexrotators.plotting as rp
from complexrotators.paths import RESULTSDIR

def main():
    PLOTDIR = os.path.join(RESULTSDIR, 'cadence_comparison')
    if not os.path.exists(PLOTDIR):
        os.mkdir(PLOTDIR)

    rp.plot_cadence_comparison(PLOTDIR, ticid='141146667', sector=48)

if __name__ == "__main__":
    main()
