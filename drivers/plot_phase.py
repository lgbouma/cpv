import os
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'phase')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

cp.plot_phase(PLOTDIR)
