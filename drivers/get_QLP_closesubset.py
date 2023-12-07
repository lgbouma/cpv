import pandas as pd, numpy as np
import os
from os.path import join

QLPDIR = "/ar1/local/QLP"
csvpath = join(QLPDIR, "QLP_s1s55_X_GDR2.csv")
df = pd.read_csv(csvpath)

PLX_CUT = 2 # mas -> 500pc
sel = (df.parallax > PLX_CUT)

sdf = df[sel]
sdf['dr2_source_id'] = sdf['dr2_source_id'].astype(str)

# get TIC IDs
qdf = pd.read_csv('/ar1/local/QLP/MyTable_v2_lukebouma_s1s55_qlp_CasJobs_TIC8.csv', dtype=str)
qdf = qdf[~pd.isnull(qdf.GAIA)]
qdf['dr2_source_id'] = qdf.GAIA.astype(str)
qdf = qdf[['ticid','dr2_source_id']]

mdf = sdf.merge(qdf, how='left', on='dr2_source_id')

csvpath = join(QLPDIR, f"QLP_s1s55_X_GDR2_parallax_gt_{PLX_CUT}.csv")
mdf.to_csv(csvpath, index=False)
print(f'Wrote {csvpath}')
