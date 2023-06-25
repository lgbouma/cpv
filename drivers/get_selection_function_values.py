import pandas as pd, numpy as np
import os
from os.path import join
from complexrotators.paths import PAPERDIR

localdir = '/Users/luke/local/SPOCLC'
csvpath = join(localdir, 'gaia_X_spoc2min_merge.csv')
df = pd.read_csv(csvpath)

sel = (
    (df.TESSMAG < 16) & (df.bp_rp > 1.5) &
    (df.M_G > 4) & (df.parallax > 1e3*(1/150))
)

sdf = df[sel]

N_1 = len(np.unique(sdf.TICID))
N_2 = len(sdf)

latex_txt = (
    r"\newcommand{\nstarssearched}{"+str(N_1)+"}\n"
    r"\newcommand{\nlcssearched}{"+str(N_2)+"}\n"
)

outpath = join(PAPERDIR, 'vals_searched_stars.tex')
with open(outpath, 'w') as f:
    f.writelines(latex_txt)
print(f"Wrote {outpath}")

