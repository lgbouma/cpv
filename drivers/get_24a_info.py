from os.path import join
import os
from glob import glob

import pandas as pd, numpy as np
from cdips.utils.catalogs import get_tic_star_information
from complexrotators.paths import RESULTSDIR, TARGETSDIR, TABLEDIR, PAPERDIR

from astropy.coordinates import SkyCoord
from astropy import units as u

from numpy import array as nparr

df = pd.read_csv(join(PAPERDIR, "table1_MRT.csv"), sep="|")

ticids = [
  "264599508", "368129164", "301676454", "141146667", "5714469", "142173958",
  "58084670", "289840926", "67897871", "442571495", "272248916", "89026133",
  "148646689", "244161191"
]

sel = df.ticid.astype(str).isin(ticids)

sdf = df[sel]

tic8_rows = []

for ix, r in sdf.iterrows():

    ticid = r['ticid']

    r_tic8 = get_tic_star_information(str(ticid))

    tic8_rows.append(r_tic8)

tic8_df = pd.concat(tic8_rows)

mdf = sdf.merge(tic8_df, left_on='ticid', right_on='ID', how='inner')

c = SkyCoord(ra=nparr(mdf['ra'])*u.deg, dec=nparr(mdf['dec'])*u.deg,
             frame='icrs')

mdf['coord'] = c.to_string('hmsdms', precision=1)

selcols = 'ticid,coord,dr3_bp_rp,period,Vmag,Jmag'.split(',')

outpath = '/Users/luke/Dropbox/Documents/proposals/2023_10_COO/CPV_24A/latex/24a_targets.tex'

mdf[selcols].to_latex(outpath, index=False, float_format="%.2f")

