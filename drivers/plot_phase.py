import os
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'phase')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

# NOTE: actual options are in the plotting script here
cp.plot_phase(PLOTDIR)
