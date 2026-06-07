import os
import complexrotators.plotting as cp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'multicolor_phase')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

cp.plot_multicolor_phase_lp12502(PLOTDIR)
cp.plot_multicolor_phase_lp12502(
    PLOTDIR,
    dates=['20231208', '20231215'],
    suffix='Dec08Dec15',
    xlim=[-0.25, 0.55],
    label_dy=1.0,
    group_dy_overrides={'20231215': [4, 4, 6, 4]},
    full_keplercam_labels=True,
)
