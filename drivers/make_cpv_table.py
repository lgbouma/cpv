"""
Construct the extended and publication versions of the 2023 catalog table.
"""

from glob import glob
from os.path import join
import pandas as pd, numpy as np
import os
from complexrotators.paths import RESULTSDIR, TABLEDIR, LITDIR
from os.path import join

from complexrotators.observability import (
    get_gaia_rows, assess_tess_holdings
)

vetdir = join(RESULTSDIR, "cpvvetter")
indir = join(TABLEDIR, "2023_catalog_table")
assert os.path.exists(indir)

def main():

    csvpath = join(indir, "20230613_LGB_RJ_uticid_quality_label.csv")
    df = pd.read_csv(csvpath, sep="|")

    rows = []
    for t in df['ticid']:
        r = get_cpvtable_row(t)
        rows.append(r)
        import IPython; IPython.embed()
        assert 0


def flatten_tdf(tdf, ticid):

    ftdf = pd.DataFrame({
        'ticid': ticid,
        'sectors': ",".join(list(tdf['sector'].astype(str))),
        'N_sectors': len(tdf),
        'N_200sec': ",".join(list(tdf['N_120sec'].astype(str))),
        'N_FFI': ",".join(list(tdf['N_FFI'].astype(str))),
    }, index=[0])

    return ftdf


def get_banyan_result(gdr2_df):

    import sys
    sys.path.append("/Users/luke/Dropbox/proj/banyan_sigma")
    from core import membership_probability

    ra, dec = float(gdr2_df.ra), float(gdr2_df.dec)
    pmra, pmdec = float(gdr2_df.pmra), float(gdr2_df.pmdec)
    epmra, epmdec = float(gdr2_df.pmra_error), float(gdr2_df.pmdec_error)
    plx, eplx = float(gdr2_df.parallax), float(gdr2_df.parallax_error)
    output = membership_probability(ra=ra, dec=dec, pmra=pmra, pmdec=pmdec,
                                    epmra=epmra, epmdec=epmdec, plx=plx, eplx=eplx,
                                    use_plx=True, use_rv=False)

    probs = np.array(output['ALL'].iloc[0].round(5))
    assocs = np.array(output['ALL'].iloc[0].index)
    banyan_df = pd.DataFrame({"prob": probs}, index=assocs)
    banyan_df['assoc'] = assocs

    agepath = join(LITDIR, "Gagne2018_table1_assoc_name_ages.csv")
    agedf = pd.read_csv(agepath, comment='#')

    banyan_df = banyan_df.merge(agedf, how='left', on='assoc')

    sel = banyan_df.prob > 1e-5
    sdf = banyan_df[sel].sort_values(by='prob', ascending=False)

    assocs = ",".join(list(sdf.assoc.astype(str)))
    probs = ",".join(list(sdf.prob.astype(str)))
    ages = ",".join(list(sdf.age.astype(str)))
    age_refs = ",".join(list(sdf.age_reference.astype(str)))

    bdf = pd.DataFrame({
        'banyan_assoc': assocs,
        'banyan_prob': probs,
        'banyan_age': ages,
        'banyan_age_refs': age_refs
    }, index=[0])

    return bdf


def get_cpvtable_row(ticid):

    gdr2_df = get_gaia_rows(ticid, allcols=1)

    tdf = assess_tess_holdings(ticid, outdir=indir)

    ftdf = flatten_tdf(tdf, ticid)

    bdf = get_banyan_result(gdr2_df)

    import IPython; IPython.embed()

    row = pd.concat((ftdf, gdr2_df, bdf), axis='columns')
    assert 0



if __name__ == "__main__":
    main()
import IPython; IPython.embed()
