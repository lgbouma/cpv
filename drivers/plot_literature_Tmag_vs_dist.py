import numpy as np, pandas as pd, matplotlib.pyplot as plt
from os.path import join
import os
from complexrotators.paths import LOCALDIR, DATADIR, TABLEDIR
from aesthetic.plot import set_style, savefig
from numpy import array as nparr

def get_detection_df():
    # good, maybe, and omit impostors
    csvpath = join(
        TABLEDIR,
        '2023_catalog_table',
        '20230613_LGB_RJ_CPV_TABLE_supplemental_selfnapplied.csv'
    )
    df = pd.read_csv(csvpath, sep="|")

    df = df[(df.quality == 0) | (df.quality == 1)]

    # all dr2...
    Tmag_pred = (
            df['phot_g_mean_mag']
            - 0.00522555 * (df['phot_bp_mean_mag'] - df['phot_rp_mean_mag'])**3
            + 0.0891337 * (df['phot_bp_mean_mag'] - df['phot_rp_mean_mag'])**2
            - 0.633923 * (df['phot_bp_mean_mag'] - df['phot_rp_mean_mag'])
            + 0.0324473
    )
    df['TESSMAG_pred'] = Tmag_pred

    return df


def get_literature_cpv_df():

    indir = join(DATADIR, 'interim')
    inpath = join(indir, "Stauffer2017_k2_class1class2_supplemented.csv")
    s17_df = pd.read_csv(inpath)
    s17_df['telescope'] = 'K2'
    s17_df['provenance'] = 'Stauffer2017'

    inpath = join(indir, "Stauffer2018_k2_class1_supplemented.csv")
    s18_df = pd.read_csv(inpath)
    s18_df['telescope'] = 'K2'
    s18_df['provenance'] = 'Stauffer2018'

    inpath = join(indir, "Stauffer2021_tess_class1class2_supplemented.csv")
    s21_df = pd.read_csv(inpath)
    s21_df['telescope'] = 'TESS'
    s21_df['provenance'] = 'Stauffer2021'

    inpath = join(indir, "Gunther_2022_tess_supplemented.csv")
    g22_df = pd.read_csv(inpath)
    g22_df['telescope'] = 'TESS'
    g22_df['provenance'] = 'Zhan2019/Gunther2022'

    # latter source is not a CPV...
    p23_df = pd.DataFrame({
        'ticid': [165184400, 141306513, 65347864],
        'distance_pc': [43.2, 50.2, 39.2],
        'TESSMAG': [12.37 ,13.65, 11.41],
        'telescope': ['TESS', 'TESS', 'TESS'],
        'provenance': 'Popinchalk2023'
    }, index=[500,501,502])

    # ptfo 8-8695
    ptfo_df = pd.DataFrame({
        'ticid': 264461976,
        'distance_pc': 349.5, # 2.8614 mas DR2
        'TESSMAG': 14.03,
        'telescope': 'TESS',
        'provenance': 'Bouma2020'
    }, index=[999])

    selcols = 'ticid,distance_pc,TESSMAG,telescope,provenance'.split(',')

    df = pd.concat((
        s17_df[selcols],
        s18_df[selcols],
        s21_df[selcols],
        g22_df[selcols],
        ptfo_df[selcols],
        p23_df[selcols]
    ))

    df = df.sort_values(by='ticid')
    df = df.reset_index(drop=True)

    # CPVs also flagged by Popinchalk2023
    p23_ticids = [
        38820496,
        201789285,
        234295610,
        206544316,
        425933644,
        425937691
    ]
    sel = df.ticid.isin(p23_ticids)
    df.loc[sel, 'provenance'] += '/Popinchalk2023'

    df = df.drop_duplicates(selcols, keep='first')
    df = df.reset_index(drop=True)

    return df


def get_twomin():

    csvpath = "/Users/luke/local/SPOCLC/gaia_X_spoc2min_merge.csv"
    tdf = pd.read_csv(csvpath)

    sel = (
        (tdf.bp_rp > 1.5) &
        (tdf.TESSMAG < 16) &
        (tdf.parallax > 1e3*(1/150)) &
        (tdf.M_G > 4) &
        (tdf.SECTOR <= 55)
    )
    # S1-S55, sel fn applied
    stdf = tdf[sel]

    return stdf



def plot_tessmag_vs_distance(showliterature=1, show120seccands=1):

    # get & clean data
    df = get_detection_df()
    df = df.sort_values(by='ticid')
    df = df.reset_index(drop=True)
    df['TESSMAG'] = df['tic8_Tmag']

    ldf = get_literature_cpv_df()
    ldf = ldf.reset_index(drop=True)

    # TODO: how many had 2 min data???
    tdf = get_twomin()
    ldf['wasInSelFn'] = ldf.ticid.isin(tdf.ID)


    # the "real detection" CPV TESSMAGS do not match the "literature CPV"
    # TESSMAGS -- they are about 0.15 to 0.2 mags fainter.  this is because the
    # "real detection" TESSMAGS come from the stassun2018 relation, while the
    # literature CPV TESSMAGS come from directly querying TIC8.

    # fix the mean
    #df['TESSMAG'] += -0.2

    # adopt the literature TESSMAGs for overlapping cases
    sel0 = df.ticid.isin(ldf.ticid)
    sel1 = ldf.ticid.isin(df.ticid)
    np.testing.assert_array_equal(
        nparr(df.loc[sel0, 'ticid']), nparr(ldf.loc[sel1, 'ticid'])
    )

    df.loc[sel0, 'TESSMAG'] = nparr(ldf.loc[sel1, 'TESSMAG'])
    df.loc[sel0, 'dist_pc'] = nparr(ldf.loc[sel1, 'distance_pc'])

    #
    # make the plot!
    #
    plt.close("all")
    set_style("science")
    fig, ax = plt.subplots(figsize=(2,2))

    # 120 second CPV detections
    sel = (df.quality == 1)
    ax.scatter(
        df[sel].dist_pc, df[sel].TESSMAG,
        c='C0', s=6, zorder=5, label='CQVs B23', linewidths=0
    )

    # 120second candidates
    sel = (df.quality == 0)
    if show120seccands:
        ax.scatter(
            df[sel].dist_pc, df[sel].TESSMAG,
            c='C0', s=6, zorder=4, label='Candidates B23',
            alpha=0.5, linewidths=0
        )


    # literature comparison
    if showliterature:
        sel = ldf.telescope == 'K2'
        ax.scatter(
            ldf[sel].distance_pc, ldf[sel].TESSMAG,
            c='gray', marker='D', s=5, linewidth=0.1, zorder=4,
            label='K2 literature'
        )

        sel = ldf.telescope == 'TESS'
        ax.scatter(
            ldf[sel].distance_pc, ldf[sel].TESSMAG,
            #c='lightgray', marker='s', s=5, linewidth=0.1, zorder=6,
            label='TESS literature',
            c='white', s=15, zorder=2, marker='o',
            #edgecolors='lightgray',
            edgecolors='k',
            linewidth=0.2, alpha=1
        )

    legend = ax.legend(loc='lower center', fontsize=5,
                       handletextpad=-0.3,
                       borderaxespad=0.8, borderpad=0.15, fancybox=True,
                       framealpha=0.8, frameon=True, ncol=2, columnspacing=0.2)
    legend.get_frame().set_linewidth(0.5)

    ax.update({
        'xlabel': 'Distance [pc]',
        'ylabel': 'TESS mag',
        'ylim': [7, 16.5],
        'xlim': [0, 200]
    })

    s = ''
    if showliterature:
        s += '_showlit'

    savefig(fig, f"../results/literaturecomp/gmag_vs_distance_CPVs{s}.png")


if __name__ == "__main__":

    plot_tessmag_vs_distance(showliterature=1)
