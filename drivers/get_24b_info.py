"""
Make an observing relevant table, given a list of TIC ID's for whatever the
upcoming semester may be.
"""
from glob import glob
from os.path import join
import pandas as pd, numpy as np, matplotlib.pyplot as plt
import os, time
from complexrotators.paths import (
    DATADIR, RESULTSDIR, TABLEDIR, LITDIR, LOCALDIR, PAPERDIR
)
from os.path import join

from complexrotators.observability import (
    get_gaia_dr2_rows, assess_tess_holdings, get_gaia_dr3_rows
)
from complexrotators.getters import get_tic8_row

from make_cpv_table import flatten_tdf

indir = join(TABLEDIR, "lit_compilation")
assert os.path.exists(indir)

def get_obstable_row(ticid, overwrite=0):

    cachecsv = join(indir, f"TIC{ticid}_cpvtable_row.csv")
    if os.path.exists(cachecsv) and not overwrite:
        return pd.read_csv(cachecsv, sep="|")

    # NOTE: a few DR2 source_id's missing->dont pull gaia data, for now.
    #gdr2_df = get_gaia_dr2_rows(ticid, allcols=1)
    #gdr3_df = get_gaia_dr3_rows(ticid)

    t8_df = get_tic8_row(ticid, indir)
    tdf = assess_tess_holdings(ticid, outdir=indir)
    ftdf = flatten_tdf(tdf, ticid)

    pd.options.display.max_rows = 5000

    row = pd.concat((ftdf, t8_df), axis='columns')

    row.to_csv(cachecsv, index=False, sep="|")
    print(f"Wrote {cachecsv}")

    return row

def main(overwrite=0):

    csvpath = join(
        indir,
        '20240304_CPV_lit_compilation_R16_S17_S18_B20_S21_Z19_G22_P23_B24.csv'
    )
    _df = pd.read_csv(csvpath, sep="|")

    rows = []
    for t in _df['ticid']:
        r = get_obstable_row(t, overwrite=overwrite)
        rows.append(r)

    rdf = pd.concat(rows).reset_index(drop=True)

    # note EPIC 204060981 / TIC 49072162 has two periods b/c two CPVs
    df = _df.merge(rdf, how='inner', on=['ticid','period_hr'])
    df = df.reset_index(drop=True)

    import IPython; IPython.embed()

if __name__ == "__main__":
    main()
