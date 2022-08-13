import os
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'phase')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

# cp.plot_phase(
#     PLOTDIR,
#     ticid='TIC_188109809',
#     lc_cadences='2min_20sec',
#     manual_period=None,
#     ylim=None,
#     binsize_minutes=10
# )

# NOTE: actual options are in the plotting script here
ticids = [
    "177309964",
    "177309964",
    "425933644",
    "425933644",
    "206544316",
    "224283342",
    "146539195"
]

for ticid in ticids:

    outdir = os.path.join(PLOTDIR, f'TIC_{ticid}')
    if not os.path.exists(outdir): os.mkdir(outdir)

    cp.plot_phase(
        outdir,
        ticid=f'TIC_{ticid}',
        lc_cadences='2min',
        ylim=None,
        binsize_minutes=10,
        t0='binmin'
    )
