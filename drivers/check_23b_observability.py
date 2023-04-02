"""
Get observing diagnostics (brightness, period, dip depth).

NOTE: to use this, you need to update the "cachedir" with logfiles under
observability.py
"""
import os
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from glob import glob
from complexrotators.observability import (
    get_gaia_rows, get_tess_stats, check_astroplan_months_observable,
    check_tesspoint, merge_to_observability_table
)
from complexrotators.paths import RESULTSDIR, TARGETSDIR, TABLEDIR
from cdips.utils.catalogs import get_tic_star_information

# begin options
listname = "20220115_RahulJayararaman_CR_list.csv"
listname = "LAH_TESS_GI_Sector_43_thru_45_CRs.csv"
listname = 'LGB_complex_rotators_ddt_merge_tic8.csv'
listname = '20230327_good_CPV_ticids_d_lt_70pc.csv'
listname = '20230401_good_CPV_ticids_d_lt_95pc.csv' # complete to 80pc, about 66% in the 80-95pc shell

site = "Palomar" # or "keck"
semester = '23b'
# end options

listid = listname.replace('.csv', '')
targetlist = os.path.join(TARGETSDIR, listname)
df = pd.read_csv(targetlist)

outdir = os.path.join(TABLEDIR, f'{semester}_observability_{listid}')
if not os.path.exists(outdir):
    os.mkdir(outdir)

for ticid in df.ticid.astype(str):

    csvpath = os.path.join(outdir, f'TIC_{ticid}.csv')

    if not os.path.exists(csvpath):
        print(42*'-')
        print(f'Beginning TIC {ticid}')
        r_gaia = get_gaia_rows(ticid)
        r_tess = get_tess_stats(ticid)
        r_tic8 = get_tic_star_information(ticid)
        r_tesspoint = check_tesspoint(
            float(r_tic8['ra']), float(r_tic8['dec']), ticid
        )
        r_astroplan = check_astroplan_months_observable(
            f'TIC_{ticid}', float(r_tic8['ra']), float(r_tic8['dec']),
            site=site
        )

        #
        # is object visible from Keck/Palomar for at least two months during the
        # semester?
        #
        # 22B: 2022/08/01 through... 2023/01/31
        # 23B: 2023/08/01 through... 2024/01/31
        okmonths = np.array([8,9,10,11,12,1])
        monthsobservable = np.array(
            r_astroplan['bestmonths'][0].split(',')
        ).astype(int)

        is_visible = len(np.intersect1d(okmonths, monthsobservable)) >= 2

        #
        # is TESS looking at the object?
        # S55: starts 08/05/22
        # S60: ends 01/18/23
        #
        # S68: start 07/29/2023
        # S74: end   01/30/2024
        #

        oktessectors = np.array([68,69,70,71,72,73,74])
        sectorsobservable = np.array(
            r_tesspoint['sector'][0].split(',')
        ).astype(int)

        is_tess_visible = len(np.intersect1d(oktessectors, sectorsobservable)) >= 1

        r_astroplan[f'is_{semester}_tess_visible'] = is_tess_visible
        r_astroplan[f'is_{semester}_{site}_visible'] = is_visible

        merge_to_observability_table(
            csvpath, r_gaia, r_tess, r_tic8, r_tesspoint, r_astroplan
        )
    else:
        print(f'Found {csvpath}')

df = pd.concat(
    (pd.read_csv(os.path.join(outdir, f'TIC_{ticid}.csv'))
     for ticid in df.ticid.astype(str))
)

df = df.reset_index()
df = df.drop(['index','ra.1','dec.1'], axis='columns')
df = df.rename({'ID':'ticid'}, axis='columns')

outpath = os.path.join(outdir, f'{semester}_CR_observability_{listid}.csv')
df.to_csv(outpath, index=False, sep='|')
print(f'Made {outpath}')

#
# make some plots
#

sdf = df[df[f"is_{semester}_{site}_visible"]]

xytuples = [
    ('ra','dec', 'linear', 'linear'),
    ('dist_pc', 'Vmag', 'linear', 'linear'),
    ('bp_rp', 'phot_g_mean_mag', 'linear', 'linear'),
    ('bp_rp', 'period', 'linear', 'linear'),
    ('bp_rp', 'period', 'linear', 'log'),
    ('a_5_95', 'period', 'log', 'log'),
]

for xy in xytuples:
    xkey, ykey = xy[0], xy[1]
    xscale, yscale = xy[2], xy[3]
    plt.close('all')
    fig, ax = plt.subplots(figsize=(4,3))
    ax.scatter(df[xkey], df[ykey], c='gray', s=2, zorder=1, label='all CRs')
    ax.scatter(sdf[xkey], sdf[ykey], c='k', s=3, zorder=2, label=f'{semester} {site} visible')
    ax.set_xlabel(xkey)
    ax.set_ylabel(ykey)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.legend(fontsize='xx-small')
    s = ''
    if xscale == 'log':
        s += '_logx'
    if yscale == 'log':
        s += '_logy'
    outpath = os.path.join(outdir, f'{xkey}_vs_{ykey}{s}.png')
    fig.savefig(outpath, bbox_inches='tight', dpi=400)
    print(f'Made {outpath}')

selcols=['ticid','bp_rp','ruwe','dist_pc','period','a_5_95','phot_g_mean_mag','Vmag','Tmag','Jmag','is_23b_tess_visible','sector']

outpath = os.path.join(outdir, f'{semester}_{site}_observable_distance_sorted_{listid}.csv')

outdf = sdf[selcols].sort_values(by='dist_pc').round(2)
outdf.to_csv(outpath, index=False, sep="|")
print(f'Made {outpath}')

import IPython; IPython.embed()
