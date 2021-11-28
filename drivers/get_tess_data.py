from os.path import join
import pandas as pd, numpy as np
from cdips_followup.quicklooktools import get_tess_data
from complexrotators.paths import DATADIR

outdir = join(DATADIR, 'photometry', 'tess')

# 2 minute targets
#df0 = pd.read_csv(join(DATADIR, 'targetlists', 'complex_rotators_ddt_merge_tic8.csv'))
df0 = pd.read_csv(join(DATADIR, 'targetlists', '20211127_EPIC_K2_LAH_objects.csv'))

for ix,r in df0.iterrows():

    try:
        ticid = str(np.int64(r['tic']))
    except ValueError:
        print(f'{ix} Skipping because no ticid...')
        continue

    if len(ticid) > 0:
        print(ix,ticid)
        _ = get_tess_data(ticid, outdir=outdir, spoc=1)

# # 20 second targets
# df1 = pd.read_csv('/Users/luke/Dropbox/Documents/proposals/2020_10_TESS_DDT/complex-rotators-ddt-20s-subset.csv')


