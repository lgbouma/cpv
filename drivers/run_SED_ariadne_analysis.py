# env: (py38_ariadne)

from os.path import join
import pandas as pd, numpy as np

from complexrotators.paths import TABLEDIR
from complexrotators.sedfit import run_SED_analysis

UNIFORMPRIORS = 1 # true if you want to impose the uniform prior set

MANUAL_TRIM_DICT = {
    # ticid: each entry in form (>xmin, <xmax, >ymin, <ymax)
    # e.g. to exclude everything about 1e-9 erg/cm2/s:
    '2234692': [ (None, None, 1e-9, None), ],
    '427460117': [ (None, None, 1e-10, None), ],
    '234295610': [ (None, None, 1e-9, None), ],
}

csvpath = join(TABLEDIR, "2023_catalog_table",
               "20230613_LGB_RJ_CPV_TABLE_supplemental_selfnapplied.csv")

df = pd.read_csv(csvpath, sep="|")

df = df.sort_values(by=['goodsectors', 'tic8_Tmag'])

ticids = df.ticid

bad_ticids = []

#DEBUG_CASES = []
#for ticid in DEBUG_CASES:
for ticid in ticids:

    trimlist = None
    if str(ticid) in MANUAL_TRIM_DICT:
        trimlist = MANUAL_TRIM_DICT[str(ticid)]

    try:
        run_SED_analysis(str(ticid), trimlist=trimlist,
                         uniformpriors=UNIFORMPRIORS)

    except Exception as e:
        print(f'Error! {e}')
        bad_ticids.append(str(ticid))

print(bad_ticids)
