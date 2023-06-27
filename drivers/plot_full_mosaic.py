import os
import complexrotators.plotting as rp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'full_lcmosaic')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

rp.plot_full_lcmosaic(PLOTDIR, sortby='Ngoodsectors_tic8_Tmag')
rp.plot_full_lcmosaic(PLOTDIR, sortby='tlc_mean_period')
assert 0
rp.plot_full_lcmosaic(PLOTDIR, sortby='tic8_Tmag')
rp.plot_full_lcmosaic(PLOTDIR, sortby='dist_pc')
