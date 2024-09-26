"""
Make an observing relevant table, given a list of TIC ID's for whatever the
upcoming semester may be.  Table contains RAs, DECs, mags, and TESS sectors
during which the sources will be visible.  (Assuming an updated version of
tesspoint).
"""
from glob import glob
from os.path import join
import pandas as pd, numpy as np, matplotlib.pyplot as plt
import os, time
from complexrotators.paths import (
    DATADIR, RESULTSDIR, TABLEDIR, LITDIR, LOCALDIR, PAPERDIR
)
from os.path import join

from complexrotators.observability import get_tess_obstable_row

indir = join(TABLEDIR, "lit_compilation")
assert os.path.exists(indir)

def main(overwrite=0):

    csvpath = join(
        indir,
        '20240304_CPV_lit_compilation_R16_S17_S18_B20_S21_Z19_G22_P23_B24.csv'
    )
    _df = pd.read_csv(csvpath, sep="|")

    rows = []
    for t in _df['ticid']:
        r = get_tess_obstable_row(t, indir, overwrite=overwrite)
        rows.append(r)

    rdf = pd.concat(rows).reset_index(drop=True)
    rdf = rdf.rename(columns={'ticid': 'ticid_x'})

    df = pd.concat((rdf,_df), axis=1)
    df.tic8_GAIA = df.tic8_GAIA.astype(str)
    outcsvpath = join(
        indir,
        '20240304_CPV_lit_compilation_R16_S17_S18_B20_S21_Z19_G22_P23_B24_TIC8_obs_supplemented.csv'
    )
    df.to_csv(outcsvpath, index=False)
    print(f'Wrote {outcsvpath}')

    max_sector = []
    for sstr in df['sectors'].astype(str):
        if sstr == '-1':
            max_sector.append(int(sstr))
        elif "," in sstr:
            max_sector.append(max(np.array(sstr.split(",")).astype(int)))
        else:
            max_sector.append(int(sstr))

    df['max_sector'] = max_sector

    #########################
    # 25A / Cycle7 specific #
    #########################
    df['c7observable'] = np.array(max_sector) >= 84
    #########################

    viable_sectors = range(88, 94+1)
    is_25a_observables = []
    for sstr in df['sectors'].astype(str):
        is_25a_observable = False
        for s in viable_sectors:
            if str(s) in sstr:
                is_25a_observable = True
                continue
        is_25a_observables.append(is_25a_observable)

    df['is_25a_observable'] = np.array(is_25a_observables).astype(int)

    selcols = [
        'ticid', 'sectors', 'c7observable', 'is_25a_observable', 'N_sectors', 'N_200sec', 'N_FFI',
        'tic8_ra', 'tic8_dec', 'tic8_Tmag',
        'tic8_GAIA',
        'tic8_plx',
        'tic8_Jmag',
        'tic8_Teff',
        'tic8_gaiabp',
        'tic8_gaiarp',
        'original_id', 'distance_pc',
        'cluster', 'bibcode',
        'period_hr', 'quality', 'telescope'
    ]
    sdf = df[selcols]
    # TRUNCRATED VERSION OF EVERYTHING
    outcsvpath = join(
        indir,
        '20240304_CPV_lit_compilation_R16_S17_S18_B20_S21_Z19_G22_P23_B24_TIC8_obs_truncated.csv'
    )
    sdf.to_csv(outcsvpath, index=False)
    print(f'Wrote {outcsvpath}')

    # OBSERVABLE IN 25a, SIMULTANEOUS WITH TESS
    outcsvpath = join(
        indir,
        '20240304_CPV_lit_compilation_R16_S17_S18_B20_S21_Z19_G22_P23_B24_TIC8_obs_truncated_25aOBSERVEABLE.csv'
    )
    sdf[sdf.is_25a_observable == 1].sort_values(by='tic8_Tmag').to_csv(
        outcsvpath, index=False
    )
    print(f'Wrote {outcsvpath}')

    # OBSERVABLE IN CYCLE7 BY TESS
    outcsvpath = join(
        indir,
        '20240304_CPV_lit_compilation_R16_S17_S18_B20_S21_Z19_G22_P23_B24_TIC8_obs_truncated_C7OBSERVABLE.csv'
    )
    sdf[sdf.c7observable].sort_values(by='tic8_Tmag').to_csv(
        outcsvpath, index=False
    )
    print(f'Wrote {outcsvpath}')

    import IPython; IPython.embed()

if __name__ == "__main__":
    main()
