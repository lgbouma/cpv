import os
import complexrotators.plotting as rp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'spectrum_windows')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

#starid = 'TIC146539195'

starids = 'TIC402980664 TIC141146667 TIC146539195 TIC264599508 TIC408188366 TIC59129133'
starids = starids.split(" ")

ylims = [
    [[0.0,35],[0.0,35],[0.3,7],[0.3,1.4]], # tic 4029
    None, # tic1411
    [[0.0,35],[0.0,35],[0.3,7],[0.3,1.4]], # tic 1465
    [[0.0,14],[0.0,14],[0.,3.5],[0.3,1.4]], # tic 2645
    None,
    None
]

for _ylims, starid in zip(ylims, starids):

    rp.plot_spectrum_windows(PLOTDIR, starid, inst='HIRES',
                             ylims=_ylims)
