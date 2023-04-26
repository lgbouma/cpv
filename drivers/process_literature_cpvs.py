"""
Get all known CPVs from Stauffer+2017, Stauffer+2018, Stauffer+2021, and
Zhan+2019/Gunther+2022.  Specifically, get their TICIDs.  Then write to a file.
(with original_identifier, alternative_identifier, and ticid columns)

Note: this approach omits PTFO 8-8695.  :-(
"""
import numpy as np, pandas as pd
import os
from os.path import join
from glob import glob
from astropy.table import Table
from astropy.io import fits
from astrobase.services.identifiers import gaiadr2_to_tic
from copy import deepcopy

from complexrotators.observability import (
    get_gaia_rows, get_tess_stats, check_astroplan_months_observable,
    check_tesspoint, merge_to_observability_table
)
from cdips.utils.catalogs import get_tic_star_information

from complexrotators.paths import DATADIR, LOCALDIR

def get_stauffer17():
    # Stauffer, K2
    s17t1_path = join(
        DATADIR, "literature",
        "Stauffer_2017_Usco_rhoOph_table1_reformatted.txt"
    )
    s17_df = pd.read_csv(s17t1_path)
    s17_df['bibcode'] = '2017AJ....153..152S'
    return s17_df

def get_stauffer18():
    s18t1_path = join(
        DATADIR, "literature",
        "Stauffer_2018_USco_rhoOph_Taurus_table1_reformatted.txt"
    )
    s18_df = pd.read_csv(s18t1_path)
    s18_df['bibcode'] = '2018AJ....155...63S'
    return s18_df

def get_stauffer21():
    # Stauffer, TESS
    s21t1_path = join(
        DATADIR, "literature",
        "Stauffer_2021_table1_ucl_lcc_scallops_reformatted.txt"
    )
    s21t1_df = pd.read_csv(s21t1_path)
    s21t1_df['cluster'] = 'USco'
    s21t1_df['Class'] = 1
    s21t1_df['bibcode'] = '2021AJ....161...60S/Table1'

    s21t2_path = join(
        DATADIR, "literature",
        "Stauffer_2021_ucl_lcc_class2_persistent_flux_dips_reformatted.txt"
    )
    s21t2_df = pd.read_csv(s21t2_path)
    s21t2_df['cluster'] = 'USco'
    s21t2_df['Class'] = 2
    s21t2_df['bibcode'] = '2021AJ....161...60S/Table2'

    s21_df = pd.concat((s21t1_df, s21t2_df))
    s21_df['ticid'] = s21_df['TICID']
    return s21_df

def get_gunther22():
    # Zhan/Gunther, TESS
    g22_path = join(
        DATADIR, "literature",
        "Gunther_2022_table1_reformatted.txt"
    )
    g22_df = pd.read_csv(g22_path)
    g22_df['bibcode'] = '2019ApJ...876..127Z and 2022AJ....163..144G'
    g22_df['ticid'] = g22_df['TICID']
    return g22_df

def get_k2_x_gaiadr2():
    k2_dr2_path = join(os.path.dirname(LOCALDIR), 'k2_dr2_1arcsec.fits')
    hl = fits.open(k2_dr2_path)
    k2_df = Table(hl[1].data).to_pandas()
    return k2_df

def left_merge(df0, df1, col0, col1):
    # execute a left-join ensuring the columns are cast as strings
    df0[col0] = df0[col0].astype(str)
    df1[col1] = df1[col1].astype(str)
    return df0.merge(
        df1, left_on=col0, right_on=col1, how='left'
    )

def get_stauffer17_xmatch(s17_df_gk=None):

    s17_csvpath = join(
        DATADIR, "literature",
        "Stauffer_2017_Usco_rhoOph_table1_x_EPIC_x_GAIADR2_x_TICID.txt"
    )
    if not os.path.exists(s17_csvpath):
        ticids = []
        for epic, s in zip(s17_df_gk.EPICID, s17_df_gk.dr2_source_id):
            if str(epic) == '204918279':
                # RIK-21 no Gaia DR2 match
                ticids.append("220765024")
            else:
                ticid = gaiadr2_to_tic(s)
                ticids.append(ticid)
        s17_df_gkt = deepcopy(s17_df_gk)
        s17_df_gkt['ticid'] = np.array(ticids).astype(str)
        s17_df_gkt.to_csv(s17_csvpath, index=False)
    else:
        s17_df_gkt = pd.read_csv(s17_csvpath)

    return s17_df_gkt


def get_stauffer18_xmatch(s18_df_gk):
    s18_csvpath = join(
        DATADIR, "literature",
        "Stauffer_2018_USco_rhoOph_Taurus_table1_x_EPIC_x_GAIADR2_x_TICID.txt"
    )
    if not os.path.exists(s18_csvpath):

        s18_df_gk.loc[
            s18_df_gk.EPICID.astype(str)=="204185983", "dr2_source_id"
        ] = "6050952515921144192"

        ticids = []
        for epic, s in zip(s18_df_gk.EPICID, s18_df_gk.dr2_source_id):
            if str(epic) == '204185983':
                # V2304 Oph no match not sure why
                ticids.append("203822419")
            else:
                ticid = gaiadr2_to_tic(s)
                ticids.append(ticid)
        s18_df_gkt = deepcopy(s18_df_gk)
        s18_df_gkt['ticid'] = np.array(ticids).astype(str)
        s18_df_gkt.to_csv(s18_csvpath, index=False)
    else:
        s18_df_gkt = pd.read_csv(s18_csvpath)

    return s18_df_gkt


def append_tessmag_given_df_with_ticid(df):
    tmags = []
    for ticid in df.ticid:
        print(ticid)
        r_tic8 = get_tic_star_information(str(ticid))
        tmags.append(float(r_tic8['Tmag']))
    df['TESSMAG'] = tmags
    return df


def main():

    #
    # k2
    #
    s17_df = get_stauffer17()
    s18_df = get_stauffer18()
    k2_df = get_k2_x_gaiadr2()
    k2_df['dr2_source_id'] = k2_df.source_id.astype(str)
    s17_df_gk = left_merge(s17_df, k2_df, 'EPICID', 'epic_number')
    s18_df_gk = left_merge(s18_df, k2_df, 'EPICID', 'epic_number')

    s17_df_gkt = get_stauffer17_xmatch(s17_df_gk)
    sel = (s17_df_gkt.Class != 3)
    s_s17_df_gkt = s17_df_gkt[sel]
    s18_df_gkt = get_stauffer18_xmatch(s18_df_gk)

    # tess
    s21_df = get_stauffer21()
    g22_df = get_gunther22()

    # what quantities do you need?
    # really, only tessmag... and distance...
    # (period is nice to have on-hand, but you dont need to clean it quite yet)

    # calculate distances
    s_s17_df_gkt['distance_pc'] = 1/ ( s_s17_df_gkt.parallax * 1e-3)
    s18_df_gkt['distance_pc'] = 1/ ( s18_df_gkt.parallax * 1e-3)
    s21_df['distance_pc'] = 1/ ( s21_df.plx * 1e-3)
    g22_df['distance_pc'] = g22_df['distance_pc']

    # get TESSMAGs
    s_s17_df_gkt = append_tessmag_given_df_with_ticid(s_s17_df_gkt)
    s18_df_gkt = append_tessmag_given_df_with_ticid(s18_df_gkt)
    s21_df = append_tessmag_given_df_with_ticid(s21_df)
    g22_df = append_tessmag_given_df_with_ticid(g22_df)

    outdir = join(DATADIR, 'interim')
    outpath = join(outdir, "Stauffer2017_k2_class1class2_supplemented.csv")
    s_s17_df_gkt.to_csv(outpath, index=False)
    print(f"Wrote {outpath}")

    outpath = join(outdir, "Stauffer2018_k2_class1_supplemented.csv")
    s18_df_gkt.to_csv(outpath, index=False)
    print(f"Wrote {outpath}")

    outpath = join(outdir, "Stauffer2021_tess_class1class2_supplemented.csv")
    s21_df.to_csv(outpath, index=False)
    print(f"Wrote {outpath}")

    outpath = join(outdir, "Gunther_2022_tess_supplemented.csv")
    g22_df.to_csv(outpath, index=False)
    print(f"Wrote {outpath}")


if __name__ == "__main__":
    main()
