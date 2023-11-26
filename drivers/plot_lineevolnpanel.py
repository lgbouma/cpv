import os
import complexrotators.plotting as rp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'lineevolnpanel')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

rp.plot_lineevolnpanel(PLOTDIR)

