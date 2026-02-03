import os, socket
from complexrotators import __path__
from os.path import join

DATADIR = join(os.path.dirname(__path__[0]), 'data')
PHOTDIR = join(DATADIR, 'photometry', 'ground')
LITDIR = join(DATADIR, 'literature')
CSVDIR = join(DATADIR, 'photometry', 'cdips_followup_csvs')
TARGETSDIR = join(DATADIR, 'targetlists')
RESULTSDIR = join(os.path.dirname(__path__[0]), 'results')
TABLEDIR = join(RESULTSDIR, 'tables')
PAPERDIR = join(os.path.dirname(__path__[0]), 'papers', 'paper')

# 120-second data is downloaded here
LKCACHEDIR = join(
    os.path.expanduser('~'), '.lightkurve', 'cache', 'mastDownload', 'TESS'
)
TARSCACHEDIR = join(
    os.path.expanduser('~'), '.tars_cache'
)

# system-dependent.
# this directory contain the following sub-directories:
#.
#├── CURL_SCRIPTS
#├── FILE_OPENS
#├── METADATA
#├── sector-1
#├── sector-11
#├── sector-{XX}...
# where the CURL_SCRIPTS are the bulk download scripts
# (tesscurl_sector_{XX}_lc.sh) from
# https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html
SPOCDIR = '/nfs/phtess2/ar0/TESS/SPOCLC'
#SPOCDIR = '/ar1/TESS/SPOCLC'
QLPDIR = '/ar1/TESS/QLP'

# banyan-sigma needs to be installed here
if socket.gethostname() == 'marduk.home':
    BANYANDIR = '/Users/luke/Dropbox/proj/banyan_sigma'
elif socket.gethostname() in ['wh1']:
    BANYANDIR = '/ar1/PROJ/luke/proj/banyan_sigma'
elif socket.gethostname() in ['wh2', 'wh3']:
    BANYANDIR = '/ar0/PROJ/luke/proj/banyan_sigma'
else:
    raise Warning('Did not find BANYANDIR; BANYAN-Σ import will fail.')

LOCALDIR = join(os.path.expanduser('~'), 'local', 'complexrotators')
for d in [LOCALDIR, TARSCACHEDIR]:
    if not os.path.exists(d): os.mkdir(d)
