import os
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'multicolor_phase')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

cp.plot_multicolor_phase(PLOTDIR)
