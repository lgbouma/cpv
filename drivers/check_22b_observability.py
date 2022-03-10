"""
Get observing diagnostics (brightness, period, dip depth). Assumes
"plot_22b_phases" has been run, and uses those cached products.
"""
import os
import pandas as pd, numpy as np
from glob import glob
from complexrotators.observability import (
    get_gaia_rows, get_tess_stats, check_astroplan_months_observable,
    check_tesspoint, merge_to_observability_table
)
from complexrotators.paths import RESULTSDIR, TARGETSDIR, TABLEDIR
from cdips.utils.catalogs import get_tic_star_information

targetlist = os.path.join(TARGETSDIR, '20220115_RahulJayararaman_CR_list.csv')
df = pd.read_csv(targetlist)

outdir = os.path.join(TABLEDIR, '22b_observability')
if not os.path.exists(outdir):
    os.mkdir(outdir)

for ticid in df.ticid.astype(str):

    csvpath = os.path.join(outdir, f'TIC_{ticid}')

    if not os.path.exists(csvpath):
        r_gaia = get_gaia_rows(ticid)
        r_tess = get_tess_stats(ticid)
        r_tic8 = get_tic_star_information(ticid)
        r_tesspoint = check_tesspoint(
            float(r_tic8['ra']), float(r_tic8['dec']), ticid
        )
        r_astroplan = check_astroplan_months_observable(
            f'TIC_{ticid}', float(r_tic8['ra']), float(r_tic8['dec']),
        )

        #
        # is object visible from Keck for at least two months during the
        # semester?
        #
        # 22B: 2022/08/01 through... 2023/01/31
        okmonths = np.array([8,9,10,11,12,1])
        monthsobservable = np.array(
            r_astroplan['bestmonths'][0].split(',')
        ).astype(int)

        is_keck_visible = len(np.intersect1d(okmonths, monthsobservable)) >= 2
        r_astroplan['is_22b_keck_visible'] = is_keck_visible

        #
        # is TESS looking at the object?
        # S55: starts 08/05/22
        # S60: ends 01/18/23
        #

        oktessectors = np.array([55,56,57,58,59,60])
        sectorsobservable = np.array(
            r_tesspoint['sector'][0].split(',')
        ).astype(int)

        is_tess_visible = len(np.intersect1d(oktessectors, sectorsobservable)) >= 1
        r_tesspoint['is_22b_tess_visible'] = is_tess_visible

        merge_to_observability_table(
            csvpath, r_gaia, r_tess, r_tic8, r_tesspoint, r_astroplan
        )

    else:
        print(f'Found {csvpath}')
