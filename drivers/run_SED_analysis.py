# env: (py38_ariadne)

from os.path import join
import pandas as pd, numpy as np

from complexrotators.paths import TABLEDIR
from complexrotators.sedfit import run_SED_analysis

csvpath = join(TABLEDIR, "2023_catalog_table",
               "20230613_LGB_RJ_CPV_TABLE_supplemental.csv")
df = pd.read_csv(csvpath, sep="|")

df = df.sort_values(by=['goodsectors', 'tic8_Tmag'])

for ticid in df.ticid:
    run_SED_analysis(str(ticid))
