import os
import complexrotators.plotting as rp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'beforeafter_mosaic')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

rp.plot_beforeafter_mosaic(PLOTDIR, sortby='tic8_Tmag')
assert 0
rp.plot_beforeafter_mosaic(PLOTDIR, sortby='Ngoodsectors_tic8_Tmag')
rp.plot_beforeafter_mosaic(PLOTDIR, sortby='tlc_mean_period')
rp.plot_beforeafter_mosaic(PLOTDIR, sortby='dist_pc')
