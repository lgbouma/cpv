"""
Get all known CPVs from the literature, and homogenize them into one table.
The table will be used for comparison of new discoveries against existing
literature; it will also be used for ground based prioritization.

"The literature" includes:
    Rebull+2016, Stauffer+2017, Stauffer+2018, Bouma+2020, Stauffer+2021,
    Zhan+2019/Gunther+2022, Popinchalk+2023, and Bouma+2024.

Columns are:
    'original_id,ticid,distance_pc,TESSMAG,cluster,bibcode,period_hr,quality,telescope'

The TIC identifer in this table is never null.

Additional documentation on this process is in /data/literature/README.txt
"""
import numpy as np, pandas as pd
import os
from os.path import join
from glob import glob
from astropy.table import Table
from astropy.io import fits
from astrobase.services.identifiers import gaiadr2_to_tic
from copy import deepcopy

from cdips.utils.catalogs import get_tic_star_information

from complexrotators.paths import DATADIR, LOCALDIR, PAPERDIR, TABLEDIR

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

    s21_df['period_hr'] = s21_df['period']*24

    s21_df['original_id'] = 'TIC ' + s21_df['TICID'].astype(str)

    return s21_df

def get_gunther22():
    # Zhan/Gunther, TESS
    g22_path = join(
        DATADIR, "literature",
        "Gunther_2022_table1_reformatted.txt"
    )
    g22_df = pd.read_csv(g22_path)
    g22_df = g22_df.rename(columns={'Cluster':'cluster', 'Prot_hr':'period_hr'})
    g22_df['bibcode'] = '2019ApJ...876..127Z and 2022AJ....163..144G'
    g22_df['ticid'] = g22_df['TICID']
    g22_df['original_id'] = 'TIC ' + g22_df['TICID'].astype(str)
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
        s17_df_gkt = pd.read_csv(s17_csvpath, dtype={'dr2_source_id':str})

    sel = (s17_df_gkt.Class != 3)
    s17_df_gkt = s17_df_gkt[sel]

    # Whether Stauffer+2017 reported multiple periods in K2
    s17_df_gkt['flag_multperiod'] = (s17_df_gkt.P2 != '---').astype(int)

    # Deal with duplicate cases.
    # Three duplicates w/ ambiguous EPIC->DR2 matches:
    # EPIC 204066898, 204364515 and 205046529; for which Stauffer2017 Table1 says:
    # UScoCTIO 80A, component B, and component B are the CPVs.
    # manually read off
    # -> DR2 6237069013821425408, DR2 6243201299401423616, and DR2 6245676471873078656 respectively
    keep_ids = [
        '6237069013821425408', '6243201299401423616', '6245676471873078656'
    ]

    dup_df = s17_df_gkt[s17_df_gkt.duplicated(subset='EPICID',keep=False)]
    nodup_df = s17_df_gkt[~s17_df_gkt.duplicated(subset='EPICID',keep=False)]
    sdup_df = dup_df[dup_df.dr2_source_id.isin(keep_ids)]
    sdup_df = sdup_df[~sdup_df.duplicated(subset='dr2_source_id', keep='first')]

    s17_df_gkt = pd.concat((nodup_df, sdup_df))

    # Get the correct CPV periods.
    periods_hr = []
    for P1, P2 in zip(s17_df_gkt['P1'], s17_df_gkt['P2']):
        if P2 == '---':
            periods_hr.append(float(P1)*24)
        else:
            if '*' in P1:
                periods_hr.append(float(P1.rstrip('*'))*24)
            elif '*' in P2:
                periods_hr.append(float(P2.rstrip('*'))*24)
    s17_df_gkt['period_hr'] = periods_hr

    # Correct RIK 21...
    s17_df_gkt.loc[s17_df_gkt.Name == 'RIK-21', 'dr2_source_id'] = '---'

    # match length to Table1 of Stauffer2017
    assert len(s17_df_gkt) == 23 - 4

    s17_df_gkt['bibcode'] = '2017AJ....153..152S'
    s17_df_gkt['cluster'] = 'USco/rhoOph'

    s17_df_gkt['original_id'] = 'EPIC ' + s17_df_gkt['EPICID'].astype(str)

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
        s18_df_gkt = pd.read_csv(s18_csvpath, dtype={'dr2_source_id':str})

    # Whether Stauffer+2017 reported multiple periods in K2
    s18_df_gkt['flag_multperiod'] = (s18_df_gkt.P2 != '---').astype(int)

    # Deal with duplicate cases.
    # Two duplicates w/ ambiguous EPIC->DR2 matches:
    # EPIC 204060981 and EPIC 246682490.
    # The former has TWO components that are each CPVs.
    # ->manually read off
    # -> DR2 6049979104536677632 and DR2 3392553852036311040
    drop_id = ['3392553852037339776']

    s18_df_gkt = s18_df_gkt[~s18_df_gkt.dr2_source_id.isin(drop_id)]

    # Get the correct CPV periods.
    periods_hr = []
    for P1, P2 in zip(s18_df_gkt['P1'], s18_df_gkt['P2']):
        if P2 == '---':
            periods_hr.append(float(str(P1).rstrip('*'))*24)
        else:
            if '*' in P1:
                periods_hr.append(float(P1.rstrip('*'))*24)
            elif '*' in P2:
                periods_hr.append(float(P2.rstrip('*'))*24)
    s18_df_gkt['period_hr'] = periods_hr

    # match length to Table1 of Stauffer2018 (+1 for EPIC 204060981)
    assert len(s18_df_gkt) == 12

    s18_df_gkt['bibcode'] = '2018AJ....155...63S'

    s18_df_gkt['original_id'] = 'EPIC ' + s18_df_gkt['EPICID'].astype(str)

    return s18_df_gkt


def append_tessmag_dist_given_df_with_ticid(df):
    tmags = []
    plxs = []
    dists = []
    for ticid in df.ticid:
        print(ticid)
        r_tic8 = get_tic_star_information(str(ticid))
        tmags.append(float(r_tic8['Tmag']))
        try:
            plxs.append(float(r_tic8['plx']))
            dists.append(1/(float(r_tic8['plx'])*1e-3))
        except TypeError:
            plxs.append(np.nan)
            dists.append(np.nan)

    df['TESSMAG'] = tmags
    df['plx'] = plxs
    df['dist_pc'] = dists

    return df


def append_tessmag_given_df_with_ticid(df):
    tmags = []
    for ticid in df.ticid:
        print(ticid)
        r_tic8 = get_tic_star_information(str(ticid))
        tmags.append(float(r_tic8['Tmag']))

    df['TESSMAG'] = tmags

    return df




def main():

    paperlist= (
        "Rebull2016,Stauffer2017,Stauffer2018,Bouma2020,Stauffer2021,"
        "Zhan2019,Gunther2022,Popinchalk2023,Bouma2024".split(',')
    )

    #
    # k2
    #
    s17_df = get_stauffer17()
    s18_df = get_stauffer18()
    k2_df = get_k2_x_gaiadr2() # 1arcsecond K2<->Gaia DR2 xmatch (Bedell)
    k2_df['dr2_source_id'] = k2_df.source_id.astype(str)
    s17_df_gk = left_merge(s17_df, k2_df, 'EPICID', 'epic_number')
    s18_df_gk = left_merge(s18_df, k2_df, 'EPICID', 'epic_number')

    s17_df_gkt = get_stauffer17_xmatch(s17_df_gk)
    s18_df_gkt = get_stauffer18_xmatch(s18_df_gk)

    # tess
    s21_df = get_stauffer21()
    g22_df = get_gunther22()

    # what quantities do you need?
    # really, only tessmag... and distance...
    # (period is nice to have on-hand, but you dont need to clean it quite yet)

    # calculate distances
    s17_df_gkt['distance_pc'] = 1/ ( s17_df_gkt.parallax * 1e-3)
    s18_df_gkt['distance_pc'] = 1/ ( s18_df_gkt.parallax * 1e-3)
    s21_df['distance_pc'] = 1/ ( s21_df.plx * 1e-3)
    g22_df['distance_pc'] = g22_df['distance_pc']

    # get TESSMAGs
    s17_df_gkt = append_tessmag_given_df_with_ticid(s17_df_gkt)
    s18_df_gkt = append_tessmag_given_df_with_ticid(s18_df_gkt)
    s21_df = append_tessmag_given_df_with_ticid(s21_df)
    g22_df = append_tessmag_given_df_with_ticid(g22_df)

    outdir = join(DATADIR, 'interim')
    outpath = join(outdir, "Stauffer2017_k2_class1class2_supplemented.csv")
    s17_df_gkt.to_csv(outpath, index=False)
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

    # Rebull+2016: from Figure6, the only really convincing one is EPIC
    # 211070495.  Others are suggestive, but could be spots.
    # Unfortunately the TESS 2minS42-S44 data on this source recover a periodic
    # signal at the right period, but don't show any interesting structure.
    r16_df = pd.DataFrame({
        'original_id': 'EPIC 211070495',
        'ticid': 353927317,
        'distance_pc': 129.1,
        'TESSMAG': 15.50,
        'cluster': 'PLE',
        'bibcode': "2016AJ....152..114R",
        'period_hr': 0.606*24,  # 14.54hr
        'quality': 1,
        'telescope': 'K2'
    }, index=[0])

    # Popinchalk+2023; latter source is not a CPV...
    p23_df = pd.DataFrame({
        'ticid': [165184400, 141306513, 65347864],
        'distance_pc': [43.2, 50.2, 39.2],
        'cluster': 'TucHor',
        'bibcode': '2023ApJ...945..114P',
        'TESSMAG': [12.37 ,13.65, 11.41],
        'telescope': 'TESS',
        'quality': [1, 1, 0],
        'period_hr': [0.66321*24, 0.55597*24, 0.32008*24],
    }, index=[500,501,502])
    p23_df['original_id'] = 'TIC ' + p23_df['ticid'].astype(str)

    ptfo_df = pd.DataFrame({
        'original_id': 'PTFO 8-8695',
        'ticid': 264461976,
        'distance_pc': 349.5, # 2.8614 mas DR2
        'cluster': 'Orion-OB1',
        'bibcode': '2020AJ....160...86B',
        'TESSMAG': 14.03,
        'telescope': 'TESS',
        'quality': 1,
        'period_hr': 10.76,
    }, index=[999])

    # Bouma+2024
    csvpath = join(PAPERDIR, 'table1_MRT.csv')
    b24_df = pd.read_csv(csvpath, sep="|", dtype={'binarityflag':str})
    b24_df['bibcode'] = '2024AJ....167...38B'
    b24_df = b24_df[b24_df.quality >= 0] # drop dippers
    b24_df = b24_df.rename(
        columns={'assoc':'cluster', 'period':'period_hr', 'tic8_Tmag':
                 'TESSMAG', 'dr3_dist_pc':'distance_pc'}
    )
    b24_df['original_id'] = 'TIC ' + b24_df['ticid'].astype(str)
    b24_df['telescope'] = 'TESS'

    # Quality flags:
    # As in Bouma+2024, except "quality 2" means "I have not independently
    # assessed this object" (which is the case for Stauffer's...).

    s17_df_gkt['quality'] = 2
    s17_df_gkt['telescope'] = 'K2'

    s18_df_gkt['quality'] = 2
    s18_df_gkt['telescope'] = 'K2'

    s21_df['quality'] = 2
    s21_df['telescope'] = 'TESS'

    g22_df['quality'] = 1
    g22_df['telescope'] = 'TESS'

    # get homogeneously formatted everything and merge...
    selcols = 'original_id,ticid,distance_pc,TESSMAG,cluster,bibcode,period_hr,quality,telescope'.split(',')

    mdf = pd.concat((
        s17_df_gkt[selcols],
        s18_df_gkt[selcols],
        ptfo_df[selcols],
        s21_df[selcols],
        g22_df[selcols],
        b24_df[selcols]
    ))

    csvpath = join(
        TABLEDIR,
        'lit_compilation',
        '20240304_CPV_lit_compilation_R16_S17_S18_B20_S21_Z19_G22_P23_B24.csv'
    )
    mdf.to_csv(csvpath, sep="|")
    print(f'Wrote {csvpath}')


if __name__ == "__main__":
    main()
