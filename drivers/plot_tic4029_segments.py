"""
LP 12-502 state-switch plots  (currently quick hack)
"""

import os
import numpy as np
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'tic4029_segments')
if not os.path.exists(PLOTDIR): os.mkdir(PLOTDIR)

cp.plot_tic4029_segments( PLOTDIR)
