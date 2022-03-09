import os
import pandas as pd
from glob import glob
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR, TARGETSDIR

PLOTDIR = os.path.join(RESULTSDIR, '22b_phase')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

ticids = ['427460117', '300651846']
manual_periods = [None, 0.343843]

for mp, ticid in zip(manual_periods, ticids):
    outdir = os.path.join(PLOTDIR, f'TIC_{ticid}')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    print(f'Beginning {ticid}')

    cp.plot_phase(
        outdir,
        ticid='TIC_'+ticid,
        lc_cadences='2min_20sec',
        binsize_minutes=10,
        MANUAL_PERIOD=mp
    )
