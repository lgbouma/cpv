"""
This test is an interesting case study in why
    '''
    import lightkurve as lk
    lk.search_lightcurve('TIC 89463560').download_all()
    '''
gives a pretty cryptic download error.
"""
from complexrotators.getters import get_2min_cadence_spoc_tess_lightcurve
import lightkurve as lk
import numpy as np

ticstr = 'TIC 89463560'
lcset = lk.search_lightcurve(ticstr)
lcc0 = lcset[(lcset.author=='SPOC') & (lcset.exptime.value==120)].download_all()
lclist0 = [l for l in lcc0]

lclist1 = get_2min_cadence_spoc_tess_lightcurve(ticstr)

for lc0,lc1 in zip(lclist0, lclist1):
    assert np.all(lc0.time == lc1.time)

print('LC download test passed!')
