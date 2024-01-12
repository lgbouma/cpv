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
PAPERDIR = join(os.path.dirname(__path__[0]), 'paper')

# 120-second data is downloaded here
LKCACHEDIR = join(
    os.path.expanduser('~'), '.lightkurve', 'cache', 'mastDownload', 'TESS'
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
CDIPSDIR = '/ar1/TESS/CDIPS'

# banyan-sigma needs to be installed here
BANYANDIR = '/ar1/PROJ/luke/proj/banyan_sigma'

LOCALDIR = join(os.path.expanduser('~'), 'local', 'complexrotators')
if not os.path.exists(LOCALDIR): os.mkdir(LOCALDIR)
