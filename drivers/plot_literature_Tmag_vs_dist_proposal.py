import numpy as np, pandas as pd, matplotlib.pyplot as plt
from os.path import join
import os
from complexrotators.paths import LOCALDIR, DATADIR, TABLEDIR, TARGETSDIR
from aesthetic.plot import set_style, savefig
from numpy import array as nparr

def get_qlp_c1_to_c4_info():

    txtpath = join(TARGETSDIR, "20240318_QLP_d_lt_100pc_CPV_labelled.txt")
    with open(txtpath, 'r') as f:
        content = f.read()

    pdfpaths = content.rstrip("\n").split()

    ticids = [os.path.basename(p).split("_")[0] for p in pdfpaths]
    sstr = [os.path.basename(p).split("_")[1] for p in pdfpaths]
    cadence = [os.path.basename(p).split("_")[2] for p in pdfpaths]

    df = pd.DataFrame({
        'ticid': ticids,
        'sectorstr': sstr,
        'cadence': cadence
    })

    from process_literature_cpvs import append_tessmag_dist_given_df_with_ticid
    df = append_tessmag_dist_given_df_with_ticid(df)

    fix_dists_dict = {
        '165688184': 25.2,
        '144980': 68.5,
        '118449916': 97.1,
        '144486786': 77.4,
        '149773560': 1/(18.1209e-3),
        '169265855': 38.4,
        '242423459': 1/(10.8702e-3),
        '272248916': 80.5,
        '34578030': 1/(7.8441e-3),
        '349307342': 1/(10.9322e-3),
        '360339486': 1/(10.0663e-3),
        '368129164': 18.3,
        '425937691': 43.1,
        '58929099': 1/(14.6065e-3), # very high RUWE though!
        '65347864': 39.2,
    }

    for k, v in fix_dists_dict.items():
        df.loc[df.ticid.astype(str)==k, 'dist_pc'] = v

    df['dist_pc'] = np.round(df.dist_pc, 1)

    return df


def plot_tessmag_vs_distance(showliterature=1, showprediction=1, basems=9,
                             f=1.5):

    # literature, via process_literature_cpvs.py
    ldf = pd.read_csv(
        join(TABLEDIR, "lit_compilation",
             "20240304_CPV_lit_compilation_R16_S17_S18_B20_S21_Z19_G22_P23_B24_TIC8_obs_truncated.csv")
    )
    ldf['TESSMAG'] = ldf.tic8_Tmag
    ldf['ticid'] = ldf.ticid.astype(str)

    # get QLP tic id's
    qlpdf = get_qlp_c1_to_c4_info()
    qlpdf['ticid'] = qlpdf.ticid.astype(str)
    qlpdf.to_csv(
        '../results/tables/20240319_qlp_0to100pc_concat.csv', index=False
    )

    selcols = 'ticid'
    ldf = ldf.drop_duplicates(selcols)
    df = qlpdf.drop_duplicates(selcols)
    df = df[df.dist_pc < 100]

    ldf = ldf.sort_values(by='ticid')
    df = df.sort_values(by='ticid')

    # correct distances
    sel = df.ticid.isin(ldf.ticid)
    sel1 = ldf.ticid.isin(df.ticid)
    np.testing.assert_array_equal(
        nparr(df.loc[sel, 'ticid']), nparr(ldf.loc[sel1, 'ticid'])
    )
    #df.loc[sel, 'TESSMAG'] = nparr(ldf.loc[sel1, 'TESSMAG'])
    df.loc[sel, 'dist_pc'] = nparr(ldf.loc[sel1, 'distance_pc'])


    #
    # make the plot!
    #
    plt.close("all")
    set_style("clean")
    fig, ax = plt.subplots(figsize=(f*2,f*2))

    # C1-C4 FFI detections
    sel = ~df.ticid.isin(ldf.ticid)
    ax.scatter(
        df[sel].dist_pc, df[sel].TESSMAG,
        c='C0', s=f*1.2*basems, zorder=5, label='New (FFI)',
        linewidths=0, marker='o'
    )
    # previously known but recovered
    ax.scatter(
        df[~sel].dist_pc, df[~sel].TESSMAG,
        c='C0', s=f*0.7*basems, zorder=5, linewidths=0, alpha=0.5,
        label='Recovered (FFI)',
    )

    # 120second candidates
    if showprediction:
        sgdf = pd.read_csv(
            '/Users/luke/Dropbox/Documents/proposals/2023_04_TESS_Cycle_VI/scripts/sgdf_yield_c6.csv'
        )
        ax.scatter(
            sgdf.dist_pc, sgdf.TESSMAG,
            c='C1', s=f*basems, zorder=4, label='Forecast (FFI)',
            alpha=0.5, linewidths=0, marker='o'
        )


    # literature comparison
    if showliterature:
        sel = ldf.telescope == 'K2'
        ax.scatter(
            ldf[sel].distance_pc, ldf[sel].TESSMAG,
            c='gray', marker='D', s=f*basems, linewidth=0.1, zorder=4,
            label='K2 lit', alpha=0.7
        )

        sel = ldf.telescope == 'TESS'
        ax.scatter(
            ldf[sel].distance_pc, ldf[sel].TESSMAG,
            #c='lightgray', marker='s', s=5, linewidth=0.1, zorder=6,
            label='TESS lit (2min)',
            c='white', s=f*1.2*basems, zorder=2, marker='o',
            #edgecolors='lightgray',
            edgecolors='k',
            linewidth=f*0.2, alpha=1
        )

    legend = ax.legend(loc='lower right', fontsize=4.8,
                       handletextpad=-0.3,
                       borderaxespad=0.8, borderpad=0.15, fancybox=True,
                       framealpha=0.8, frameon=True, ncol=2, columnspacing=0.2)
    legend.get_frame().set_linewidth(0.5)

    ax.update({
        'xlabel': 'Distance [pc]',
        'ylabel': '$T$ [mag]',
        'ylim': [7, 16.5],
        'xlim': [0, 200]
    })

    s = ''
    if showliterature:
        s += '_showlit'

    # manuscript default
    if basems == 9:
        savefig(fig, f"../results/literaturecomp/proposal_tmag_vs_distance_CPVs{s}.png")
    else:
        savefig(
            fig,
            f"../results/literaturecomp/proposal_tmag_vs_distance_CPVs{s}_basems{basems}.png"
        )


if __name__ == "__main__":

    # sales plots
    plot_tessmag_vs_distance(showliterature=1, basems=12, f=1.1)

    # manuscript
    plot_tessmag_vs_distance(showliterature=1)
