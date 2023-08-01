import os
import complexrotators.plotting as rp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'spectrum_windows')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

#starid = 'TIC146539195'

starids = 'TIC141146667 TIC146539195 TIC264599508 TIC408188366 TIC59129133'
starids = starids.split(" ")

for starid in starids:
    rp.plot_spectrum_windows(PLOTDIR, starid, inst='HIRES')

