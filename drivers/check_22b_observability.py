"""
Get observing diagnostics (brightness, period, dip depth). Assumes
"plot_22b_phases" has been run, and uses those cached products.
"""
import os
import pandas as pd
from glob import glob
from complexrotators.observability import (
    get_gaia_rows, get_tess_stats, check_astroplan_observability,
    check_tesspoint, merge_to_obsevability_table
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
        r_astroplan = check_astroplan_observability(ticid)
        import IPython; IPython.embed()
        assert 0
        merge_to_obsevability_table(
            csvpath, r_gaia, r_tesspoint, r_astroplan, r_tesspoint
        )
    else:
        print(f'Found {csvpath}')
