"""
Construct the extended and publication versions of the 2023 catalog table.

Uses log output created by running find_CPVs.py using
sample_id=="2023catalog_LGB_RJ_concat" first (see HOWTO.md)
"""

from glob import glob
from os.path import join
import pandas as pd, numpy as np
import os
from complexrotators.paths import RESULTSDIR, TABLEDIR, LITDIR, LOCALDIR
from os.path import join

from astroquery.mast import Catalogs

from complexrotators.observability import (
    get_gaia_rows, assess_tess_holdings
)
from complexrotators import pipeline_utils as pu

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


def get_tic8_row(t):

    ticstr = f"TIC {t}"
    tic_data = Catalogs.query_object(ticstr, catalog="TIC")

    t8_row = pd.DataFrame(tic_data.to_pandas().iloc[0]).T

    t8_row = t8_row.rename({c:f"tic8_{c}" for c in t8_row.columns},
                           axis='columns')

    return t8_row


def get_tess_lc_properties(ticid):

    sample_id = "2023catalog_LGB_RJ_concat"
    cachedir = join(LOCALDIR, "cpv_finding")
    cachedir = join(cachedir, sample_id)

    cand_logpaths = glob(join(cachedir, f"tess*00{ticid}-*runstatus.log"))
    foundexitcodes = []

    msg = f'TIC{ticid} got no logs!!'
    assert len(cand_logpaths) > 0, msg

    sectors, periods, pgconds, n_peaks, peak_phaseunits, peak_props, nspls = (
        [], [], [], [], [], [], []
    )
    for cand_logpath in cand_logpaths:
        st = pu.load_status(cand_logpath)

        stp = st['cpv_periodsearch_results']
        stc = st['count_phased_local_minima_results']

        sector = int(stp['sector'])
        period = float(stp['period'])
        periodogram_condition = bool(stp['periodogram_condition'])

        _n_peaks = int(stc['n_peaks'])
        _peaks_phaseunits = list(np.array(eval(stc['peaks_phaseunits'])).round(2))

        array = np.array
        properties = eval(stc['properties'])

        nsplines_singlephase = stc['nsplines_singlephase']

        sectors.append(sector)
        periods.append(period)
        pgconds.append(periodogram_condition)
        n_peaks.append(_n_peaks)
        peak_phaseunits.append(_peaks_phaseunits)
        peak_props.append(properties)
        nspls.append(nsplines_singlephase)

    tlc_df = pd.DataFrame({
        'tlc_mean_period': np.nanmean(periods).round(6),
        'tlc_mean_npeaks': np.nanmean(n_peaks).round(1),
        'tlc_sectors': ",".join(np.array(sectors).astype(str)),
        'tlc_npeaks': ",".join(np.array(n_peaks).astype(str)),
        'tlc_periods': ",".join(np.array(periods).round(5).astype(str)),
        'tlc_peakloc_phaseunits': repr(peak_phaseunits).replace("],", "];")[1:-1]
    }, index=[0])

    return tlc_df



def get_cpvtable_row(ticid):

    gdr2_df = get_gaia_rows(ticid, allcols=1)

    tdf = assess_tess_holdings(ticid, outdir=indir)

    ftdf = flatten_tdf(tdf, ticid)

    bdf = get_banyan_result(gdr2_df)

    t8_df = get_tic8_row(ticid)

    tlc_df = get_tess_lc_properties(ticid)

    pd.options.display.max_rows = 5000

    row = pd.concat((ftdf, gdr2_df, bdf, t8_df, tlc_df), axis='columns')

    # TODO : need to add these
    bonus_cols = (
        "teff_sedfit rstar_sedfit temp_suppression radius_inflation "
        "teff_noactivity rstar_noactivity mstar_baraffe2015 logg_baraffe2015 "
        "R_corotation t_cross s_max manual_notes R10k_spectra_available"
    ).split()
    for c in bonus_cols:
        row[c] = 'PENDING'

    print(row.T)

    return row



if __name__ == "__main__":
    main()
import IPython; IPython.embed()
