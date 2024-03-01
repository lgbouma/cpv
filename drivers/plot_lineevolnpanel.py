import os
import complexrotators.plotting as rp
from complexrotators.paths import RESULTSDIR

jstr = 'j537'
starid = 'TIC141146667'
PLOTDIR = os.path.join(RESULTSDIR, f'lineevolnpanel_{starid}_{jstr}')
if not os.path.exists(PLOTDIR): os.mkdir(PLOTDIR)
rp.plot_lineevolnpanel(PLOTDIR, starid=starid, jstr=jstr)
assert 0


jstr = 'j531'
starid = 'TIC402980664'
PLOTDIR = os.path.join(RESULTSDIR, f'lineevolnpanel_{starid}_{jstr}')
if not os.path.exists(PLOTDIR): os.mkdir(PLOTDIR)
rp.plot_lineevolnpanel(PLOTDIR, starid=starid, jstr=jstr)


