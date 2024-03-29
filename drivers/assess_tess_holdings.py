"""
Given a list of ticids, assess what TESS data are available for them.

a) make a table where each row is a sector/cam/ccd tuple, with columns of

'ticid', 'baseline_sector_max55', 'sector', 'cam', 'ccd', 'colpix', 'rowpix',
'N_20sec', 'N_120sec', 'N_200sec', 'N_600sec', 'N_1800sec', 'N_FFI'

b) trim this table to stars with an interesting baseline duration for CPV
purposes (eg. "at least two years" or similar).

extra notes at
~/Dropbox/proj/cpv/doc/20230530_tess_holdings_archival_data_assessment.txt
"""
import os
from os.path import join
from glob import glob
import numpy as np, pandas as pd, matplotlib.pyplot as plt

from complexrotators.paths import TARGETSDIR, TABLEDIR
from complexrotators.observability import assess_tess_holdings

def wrap_assess_tess_holdings(sample_id):

    csvpath = join(TARGETSDIR, f"{sample_id}.csv")

    df = pd.read_csv(csvpath)

    outdir = join(TABLEDIR, f"tess_holdings_{sample_id}")
    if not os.path.exists(outdir): os.mkdir(outdir)

    for t in df['ticid']:
        _ = assess_tess_holdings(t, outdir)

    csvpath = join(TARGETSDIR, f"{sample_id}.csv")
    df = pd.read_csv(csvpath)

    csvpaths = glob(join(outdir, 'TIC*csv'))

    rdf = pd.concat((pd.read_csv(f) for f in csvpaths))

    outcsv = join(outdir, f"summarized_tess_holdings_20230530_{sample_id}.csv")
    rdf.to_csv(outcsv, index=False)

    # stars which were on silicon in <=sector 55, but which did not have any
    # 2-minute data
    SECTORMAX = 55
    missing2min_ticids = rdf[(rdf.N_120sec == 0) & (rdf.sector <= SECTORMAX)].ticid

    missingffi_ticids = rdf[(rdf.N_FFI == 0) & (rdf.sector <= SECTORMAX)].ticid

    pd.options.display.max_rows = 5000

    N_missing = len(np.unique(missing2min_ticids))
    N_missing_ffi = len(np.unique(missingffi_ticids))
    N_total = len(np.unique(rdf.ticid))

    print(rdf[rdf.ticid.isin(missing2min_ticids)])

    print(f'Wrote {outcsv}...')
    print(f'{N_missing} stars obsd in <=S{SECTORMAX}, but DO NOT HAVE 2-MINUTE DATA')
    print(f'{N_missing_ffi} stars obsd in <=S{SECTORMAX}, and DO NOT HAVE FFI data (currently on MAST HLSPs)')
    print(f'...out of {N_total} total good+maybe CPVs')

    # real baseline counts...
    # may30, 2023 max sector is 64
    srdf = rdf[rdf.sector <= 64].reset_index()
    fn = lambda df: max(df.sector) - min(df.sector)
    baseline_sector_max64 = np.array(srdf.groupby(srdf.ticid).apply(fn))

    srdf = rdf[rdf.sector <= 55].reset_index()
    fn = lambda df: max(df.sector) - min(df.sector)
    baseline_sector_max55 = np.array(srdf.groupby(srdf.ticid).apply(fn))

    plt.close('all')
    plt.hist(baseline_sector_max64, bins=np.arange(0,68), label='max s64',
             alpha=0.5)
    plt.hist(baseline_sector_max55, bins=np.arange(0,68), label='max s55',
             alpha=0.5)
    plt.xlabel('max(sector) - min(sector)')
    plt.ylabel('count')
    plt.ylim((0,25))
    plt.legend()
    plt.title('real baselines')
    plt.savefig(join(outdir, 'hist_real_baselines.png'), dpi=300,
                bbox_inches='tight')

    # baselines requiring only 2 minute data
    srdf = rdf[(rdf.sector <= 64) & (rdf.N_120sec != 0)].reset_index()
    fn = lambda df: max(df.sector) - min(df.sector)
    baseline_sector_max64 = np.array(srdf.groupby(srdf.ticid).apply(fn))

    srdf = rdf[(rdf.sector <= 55) & (rdf.N_120sec != 0)].reset_index()
    fn = lambda df: max(df.sector) - min(df.sector)
    baseline_sector_max55 = np.array(srdf.groupby(srdf.ticid).apply(fn))

    plt.close('all')
    plt.hist(baseline_sector_max64, bins=np.arange(0,68), label='max s64',
             alpha=0.5)
    plt.hist(baseline_sector_max55, bins=np.arange(0,68), label='max s55',
             alpha=0.5)
    plt.ylim((0,25))
    plt.legend()
    plt.xlabel('max(sector) - min(sector)')
    plt.ylabel('count')
    plt.title('baselines with only 2min data')
    plt.savefig(join(outdir, 'hist_2min_baselines.png'), dpi=300,
                bbox_inches='tight')

    # get the actual stars w/ the long baselines from 2min data
    _df = pd.DataFrame({
        'ticid': srdf.groupby(srdf.ticid).apply(fn).index,
        'baseline_sector_max55': baseline_sector_max55
    })

    mdf = _df.merge(rdf, how='left', on='ticid')

    mdf_before_after = mdf[mdf['baseline_sector_max55'] >= 20]

    outcsv = join(outdir, f"before_after_stars_S55_max_2mindataonly_20230530_{sample_id}.csv")
    mdf_before_after.to_csv(outcsv, index=False)
    print(f"Wrote {outcsv}")

    # get the subset that are known "good CPVs"
    csvpath = join(TARGETSDIR, f"20230411_good_CPV_ticids_d_lt_150pc.csv")
    gdf = pd.read_csv(csvpath)

    is_good_CPV = mdf_before_after.ticid.isin(gdf.ticid)
    has_2min_today = mdf_before_after.N_120sec > 0

    outcsv = join(outdir, f"before_after_stars_S55_max_2mindataonly_goodCPVonly_20230530_{sample_id}.csv")
    mdf_before_after[is_good_CPV & has_2min_today].to_csv(outcsv, index=False)
    print(f"Wrote {outcsv}")


if __name__ == "__main__":
    sample_id = '20230411_goodandmaybe_CPV_ticids_d_lt_150pc'
    wrap_assess_tess_holdings(sample_id)
