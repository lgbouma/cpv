import os, socket
from complexrotators import __path__

DATADIR = os.path.join(os.path.dirname(__path__[0]), 'data')
PHOTDIR = os.path.join(DATADIR, 'photometry', 'ground')
CSVDIR = os.path.join(DATADIR, 'photometry', 'cdips_followup_csvs')
TARGETSDIR = os.path.join(DATADIR, 'targetlists')
RESULTSDIR = os.path.join(os.path.dirname(__path__[0]), 'results')
TABLEDIR = os.path.join(RESULTSDIR, 'tables')

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
SPOCDIR = '/ar1/TESS/SPOCLC'

LOCALDIR = os.path.join(os.path.expanduser('~'), 'local', 'complexrotators')
if not os.path.exists(LOCALDIR): os.mkdir(LOCALDIR)
