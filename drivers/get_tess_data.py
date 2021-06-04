import os
import pandas as pd, numpy as np
from cdips_followup.quicklooktools import get_tess_data
from complexrotators.paths import DATADIR

outdir = os.path.join(DATADIR, 'photometry', 'tess')

# 2 minute targets
df0 = pd.read_csv('/Users/luke/Dropbox/Documents/proposals/2020_10_TESS_DDT/complex_rotators_ddt_merge_tic8.csv')

for ix,r in df0.iterrows():

    ticid = str(r['tic'])
    print(ix,ticid)
    _ = get_tess_data(ticid, outdir=outdir, spoc=1)

# # 20 second targets
# df1 = pd.read_csv('/Users/luke/Dropbox/Documents/proposals/2020_10_TESS_DDT/complex-rotators-ddt-20s-subset.csv')


