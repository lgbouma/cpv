import os, socket
from complexrotators import __path__

DATADIR = os.path.join(os.path.dirname(__path__[0]), 'data')
PHOTDIR = os.path.join(DATADIR, 'photometry', 'ground')
RESULTSDIR = os.path.join(os.path.dirname(__path__[0]), 'results')

LOCALDIR = os.path.join(os.path.expanduser('~'), 'local', 'complexrotators')
if not os.path.exists(LOCALDIR):
    os.mkdir(LOCALDIR)
