import os
import complexrotators.plotting as rp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'spectrum_windows')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

rp.plot_spectrum_windows(PLOTDIR)

