import pandas as pd, numpy as np
import os
from os.path import join
from complexrotators.paths import PAPERDIR

def add_comma(num):
	# latexify numbers greater than 1{,}000
    num_str = str(num)
    if num < 1000:
        return num_str
    else:
        return add_comma(num // 1000) + '{,}' + num_str[-3:]

localdir = '/Users/luke/local/SPOCLC'
csvpath = join(localdir, 'gaia_X_spoc2min_merge.csv')
df = pd.read_csv(csvpath)

sel = (
    (df.TESSMAG < 16) & (df.bp_rp > 1.5) &
    (df.M_G > 4) & (df.parallax > 1e3*(1/150)) &
    (df.SECTOR <= 55)
)

sdf = df[sel]

N_1 = len(np.unique(sdf.TICID))
N_2 = len(sdf)

latex_txt = (
    r"\newcommand{\nstarssearched}{"+add_comma(N_1)+"}\n"
    r"\newcommand{\nlcssearched}{"+add_comma(N_2)+"}\n"
)

outpath = join(PAPERDIR, 'vals_searched_stars.tex')
with open(outpath, 'w') as f:
    f.writelines(latex_txt)
print(f"Wrote {outpath}")

