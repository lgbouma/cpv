import os
import complexrotators.plotting as rp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'lc_mosaic')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

rp.plot_lc_mosaic(PLOTDIR)

#TODO FIXME FINISH THE MOSAIC in an evening or sometime when you need two brain cells
