import os
import complexrotators.plotting as rp
from complexrotators.paths import RESULTSDIR

PLOTDIR = os.path.join(RESULTSDIR, 'winered_windows')
if not os.path.exists(PLOTDIR):
    os.mkdir(PLOTDIR)

frameids = [
    'WINA00038979',
    'WINA00038982',
    'WINA00038983',
    'WINA00038984'
]
for frameid in frameids:
    rp.plot_winered_windows(PLOTDIR, inst='WINERED', datestr='june10',
                            frameid=f'{frameid}', starid='TIC_167664935')

for ix in range(475, 484, 1):
    rp.plot_winered_windows(PLOTDIR, inst='WINERED', datestr='june03',
                            frameid=f'WINA00037{ix}', starid='TIC_89026133')

for ix in range(986, 993, 1):
    rp.plot_winered_windows(PLOTDIR, inst='WINERED', datestr='june10',
                            frameid=f'WINA00038{ix}', starid='TIC_89026133')



