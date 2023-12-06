"""
Step9&11 of the HOWTO.
"""
import pandas as pd, numpy as np
import os
from os.path import join
from cdips.utils.gaiaqueries import given_votable_get_df

QLPDIR = "/Users/luke/local/QLP"

gzpath = join(QLPDIR, "qlp_s1s55_firsthalf-result.vot.gz")
gzpath2 = join(QLPDIR, "qlp_s1s55_firsthalf-result.vot.gz")

if not os.path.exists(gzpath):

    csvpath = os.path.join(
        QLPDIR, "MyTable_v2_lukebouma_s1s55_qlp_CasJobs_TIC8.csv"
    )

    df = pd.read_csv(csvpath, dtype=str)

    # some subset of TIC doesnt have Gaia Xmatches.  Drop them.
    sel = ~pd.isnull(df.GAIA)
    sdf = df[sel]

    outpath = os.path.join(
        QLPDIR, "MyTable_v2_lukebouma_s1s55_qlp_CasJobs_TIC8_GaiaDR2ids.csv"
    )
    sdf['GAIA'].to_csv(outpath, index=False)

    print("Gotta run the Gaia archive query")
    assert 0

gdf1 = given_votable_get_df(gzpath)
gdf2 = given_votable_get_df(gzpath2)

gdf = pd.concat((gdf1, gdf2))

gdf['M_G'] = gdf['phot_g_mean_mag'] + 5*np.log10(gdf['parallax']/1e3) + 5

outpath = join(QLPDIR, "QLP_s1s55_X_GDR2.csv")
gdf.to_csv(outpath, index=False)
