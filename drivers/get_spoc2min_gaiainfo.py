"""
Step9&11 of the HOWTO.
"""
import pandas as pd, numpy as np
import os
from os.path import join
from cdips.utils.gaiaqueries import given_votable_get_df

SPOCDIR = "/Users/luke/local/SPOCLC"

gzpath = join(SPOCDIR, "s1s58_sample-result.vot.gz")

if not os.path.exists(gzpath):

    csvpath = os.path.join(SPOCDIR, "MyTable_lukebouma_s1_s58_CasJobs_TIC8.csv")

    df = pd.read_csv(csvpath, dtype=str)

    # N=445444 unique TICIDs with 2-minute data in S1-S58
    # N=2440 of them do not have Gaia IDs.  (I wonder why?)
    sel = ~pd.isnull(df.GAIA)
    sdf = df[sel]

    outpath = os.path.join(
        SPOCDIR, "MyTable_lukebouma_s1_s58_CasJobs_TIC8_GAIADR2_ids.csv"
    )
    sdf['GAIA'].to_csv(outpath, index=False)

    print("Gotta run the Gaia archive query")

gdf = given_votable_get_df(gzpath)

gdf['M_G'] = gdf['phot_g_mean_mag'] + 5*np.log10(gdf['parallax']/1e3) + 5

outpath = join(SPOCDIR, "spoc2min_s1_s58_X_GDR2.csv")
gdf.to_csv(outpath, index=False)
