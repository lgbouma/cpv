"""
Step 12 of the HOWTO process of getting all the SPOC 2-minute LC metadata.

Merges all the CSV files into one large one, appropriate for messing around in
glue.
"""

import pandas as pd, numpy as np
from glob import glob
from os.path import join
import os

SPOCDIR = '/ar1/TESS/SPOCLC'

csvpaths = glob(join(SPOCDIR, "METADATA", "spoc2min_*_sector-*.csv"))

df = pd.concat((pd.read_csv(f) for f in csvpaths))

csvpaths = glob(join(SPOCDIR, "METADATA", "spoc2min_*_sector-*.csv"))

outpath = join(SPOCDIR, "METADATA", "MERGED_spoc2min_sector1_to_sector58.csv")

df.to_csv(outpath, index=False)
print(f"Wrote {outpath}")
