"""
Given a TIC ID, make a river plot.
"""

from complexrotators.paths import RESULTSDIR

ticid = '177309964'

from astrobase.services.tesslightcurves import (
    get_two_minute_spoc_lightcurves
)

#FIXME

outdir = ''
lcfiles = glob(os.path.join(outdir,'mastDownload','TESS','*','tess*fits'))
if len(lcfiles) == 0:
    lcfiles = get_two_minute_spoc_lightcurves(ticid, download_dir=outdir)

