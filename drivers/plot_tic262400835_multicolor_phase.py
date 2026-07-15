import os
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'multicolor_phase')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

cp.plot_multicolor_phase(PLOTDIR)

# variant using the 2026-07 Murgas MuSCAT2 re-extraction (.dat light curves);
# saves TIC262400835_multicolor_phase_stacked_reextracted.png. Drop the first
# and last 30 min of the 2020 Dec 15 night to remove its edge-of-night offset.
cp.plot_multicolor_phase(PLOTDIR, reextracted=True, trim_edges_min=30)
