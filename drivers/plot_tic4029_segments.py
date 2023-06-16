"""
LP 12-502 state-switch plots  (currently quick hack)
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from astrobase.lcmath import (
    phase_magseries, phase_bin_magseries, sigclip_magseries,
    find_lc_timegroups, phase_magseries_with_errs, time_bin_magseries
)
from aesthetic.plot import set_style


# sector 25 showed some 
#hl = fits.open("/Users/luke/.lightkurve-cache/mastDownload/TESS/tess2020133194932-s0025-0000000402980664-0182-s/tess2020133194932-s0025-0000000402980664-0182-s_lc.fits")
# sector 26
hl = fits.open("/Users/luke/.lightkurve-cache/mastDownload/TESS/tess2020160202036-s0026-0000000402980664-0188-s/tess2020160202036-s0026-0000000402980664-0188-s_lc.fits")


d = hl[1].data

hl.close()

bd = time_bin_magseries(d['TIME'], d['PDCSAP_FLUX'], binsize=1200, minbinelems=1)


set_style("science")
fig, ax = plt.subplots(figsize=(12,4))
ax.plot(bd['binnedtimes'], bd['binnedmags'], c='k', zorder=2)

plt.show()
