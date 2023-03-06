import pandas as pd, numpy as np
import os
from os.path import join
from cdips.utils.gaiaqueries import given_votable_get_df

def left_merge(df0, df1, col0, col1):
    # execute a left-join ensuring the columns are cast as strings
    df0[col0] = df0[col0].astype(str)
    df1[col1] = df1[col1].astype(str)
    return df0.merge(
        df1, left_on=col0, right_on=col1, how='left'
    )

SPOCDIR = "/nfs/phtess2/ar0/TESS/SPOCLC"

# gaiadr2 source_id, ticid,
df = pd.read_csv(
    join(SPOCDIR, "MyTable_lukebouma_s1_s58_CasJobs_TIC8.csv"),
    dtype=str
)
# TESS 2min header info
df_t = pd.read_csv(
    join(SPOCDIR, "MERGED_spoc2min_sector1_to_sector58.csv")
)
# Gaia info
df_g = pd.read_csv(
    join(SPOCDIR, "spoc2min_s1_s58_X_GDR2.csv")
)

mdf0 = left_merge(df_g, df, "dr2_source_id", "GAIA")

mdf = left_merge(mdf0, df_t, "uticid", "TICID")

dupcols = ["TICID", "SECTOR", "CAMERA", "CCD", "fitspath"]

# a few oddball cases (0.2%)
mdf = mdf.drop_duplicates(dupcols)

outcsv = join(SPOCDIR, "gaia_X_spoc2min_merge.csv")
mdf.to_csv(outcsv, index=False)
print(f"Wrote {outcsv}")
