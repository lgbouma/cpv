"""
Collect manual labelling results from /results/cpvvetter into the list of CPVs
to be used in a catalog paper.
"""

from glob import glob
from os.path import join
import pandas as pd, numpy as np
import os
from complexrotators.paths import RESULTSDIR
from os.path import join

vetdir = join(RESULTSDIR, "cpvvetter")

####################################
# get labels from rahul's pipeline #
####################################
r_fpaths = glob(join(vetdir, "rahul_20230501_labelled", "*pdf"))
r_fnames = [os.path.basename(f) for f in r_fpaths]
ticids = [f.split("_")[0] for f in r_fnames]
sectorstr = [f.split("_")[1] for f in r_fnames]
labels = [f.split("[")[-1].split("]")[0] for f in r_fnames]

r_df = pd.DataFrame({
    'ticid': ticids,
    'sectorstr': sectorstr,
    'label': labels,
    'pipeline': 'Fourier',
    'fpath': r_fpaths,
})

##################################
# get labels from lgb's pipeline #
##################################
l_fpaths = glob(join(vetdir, "*mkdwarf*labelled*", "*pdf"))
l_fnames = [os.path.basename(f) for f in l_fpaths]
ticids = [f.split("_")[0] for f in l_fnames]
sectorstr = [f.split("_")[1] for f in l_fnames]
labels = [f.split("[")[-1].split("]")[0] for f in l_fnames]

l_df = pd.DataFrame({
    'ticid': ticids,
    'sectorstr': sectorstr,
    'label': labels,
    'pipeline': 'DipCount',
    'fpath': l_fpaths,
})

#############################################################
# merge.  anything labelled "good" in any sector is "good". #
# ignore anything that isn't either a "good" or "maybe" CPV.#
#############################################################

mdf = pd.concat((r_df, l_df))

mdf = mdf.sort_values(by=['ticid', 'sectorstr'])

sel_CPV = mdf['label'].str.contains("CPV")

mdf = mdf[sel_CPV]

# no "ratty" objects
sel = mdf['label'].str.contains("ratty")
assert len(mdf[sel]) == 0

# CPV quality status is either "good" or "maybe", never both.
sel = (mdf['label'].str.contains("good")) & (mdf['label'].str.contains("maybe"))
assert len(mdf[sel]) == 0

sel = (mdf['label'].str.contains("good")) | (mdf['label'].str.contains("maybe"))
assert len(mdf[sel]) == len(mdf)

# sort for good first, maybe second
_df0 = mdf[mdf["label"].str.contains("good")]
_df1 = mdf[mdf["label"].str.contains("maybe")]
assert len(pd.concat((_df0, _df1))) == len(mdf)
mdf = pd.concat((_df0, _df1))

from complexrotators.paths import TABLEDIR
outdir = join(TABLEDIR, "2023_catalog_table")
if not os.path.exists(outdir): os.mkdir(outdir)
csvpath = join(outdir, "20230613_LGB_RJ_good_maybe_concat_sectorlevel.csv")
mdf.to_csv(csvpath, index=False)
print(f"Wrote {csvpath}")

# now assemble a dataframe where each row is a unique TIC ID, and there are
# lists of "good" sectors, vs "maybe" sectors
u_ticids = np.unique(mdf.ticid)

goodsectors, maybesectors = [], []

for ticid in u_ticids:

    r = mdf[mdf.ticid == ticid]

    s_good = r['label'].str.contains("good")
    s_maybe = r['label'].str.contains("maybe")

    good_sectors = np.array(
        [s[1:].lstrip("0") for s in np.unique(r.loc[s_good, 'sectorstr'])]
    )
    maybe_sectors = np.array(
        [s[1:].lstrip("0") for s in np.unique(r.loc[s_maybe, 'sectorstr'])]
    )
    maybe_sectors = maybe_sectors[~np.in1d(maybe_sectors, good_sectors)]

    goodsectors.append(",".join(good_sectors))
    maybesectors.append(",".join(maybe_sectors))

outdf = pd.DataFrame({
    'ticid': u_ticids,
    'goodsectors': goodsectors,
    'maybesectors': maybesectors
})
outdf['goodsectors'] = outdf.goodsectors.astype(str)
outdf['maybesectors'] = outdf.maybesectors.astype(str)
csvpath = join(outdir, "20230613_LGB_RJ_uticid_quality_label.csv")
outdf.to_csv(csvpath, index=False, sep="|")
print(f"Wrote {csvpath}")
