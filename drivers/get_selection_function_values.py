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

from getters import get_cqv_search_sample
sdf = get_cqv_search_sample()

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

