"""
Construct the extended and publication versions of the 2023 catalog table.

Uses log output created by running find_CPVs.py using
sample_id=="2023catalog_LGB_RJ_concat" first (see HOWTO.md)

Utilities:

| get_cpvtable_row(ticid)
| gdr2_df = get_gaia_rows(ticid)
| tdf = assess_tess_holdings(ticid)
| ftdf = flatten_tdf(tdf, ticid)
| bdf = get_banyan_result(gdr2_df)
| t8_df = get_tic8_row(ticid)
| tlc_df = get_tess_cpv_lc_properties(ticid)
| sed_df = get_sedfit_results(ticid)
| iso_df = get_isochrone_mass(sed_df, bdf)

"""

from glob import glob
from os.path import join
import pandas as pd, numpy as np, matplotlib.pyplot as plt
import os, time
from complexrotators.paths import (
    DATADIR, RESULTSDIR, TABLEDIR, LITDIR, LOCALDIR, PAPERDIR
)
from os.path import join

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from astroquery.mast import Catalogs
from astroquery.exceptions import ResolverError

from complexrotators.observability import (
    get_gaia_rows, assess_tess_holdings
)
from complexrotators.isochroneinterp import get_PARSEC
from complexrotators import pipeline_utils as pu

vetdir = join(RESULTSDIR, "cpvvetter")
indir = join(TABLEDIR, "2023_catalog_table")
assert os.path.exists(indir)

def main(overwrite=0):

    csvpath = join(indir, "20230613_LGB_RJ_uticid_quality_label.csv")
    _df = pd.read_csv(csvpath, sep="|")

    csvpath2 = join(DATADIR, "targetlists", "20230813_DEBUNKED_CQVs.txt")
    debunked_df = pd.read_csv(csvpath2)

    csvpath3 = join(DATADIR, "targetlists", "20230411_goodandmaybe_CPV_ticids_d_lt_150pc.csv")
    dipcount_df = pd.read_csv(csvpath3)

    csvpath4 = join(DATADIR, "targetlists", "20230501_RAHUL_FULL_LIST_NO_DUPLICATES.csv")
    fourier_df = pd.read_csv(csvpath4)

    rows = []
    for t in _df['ticid']:
        r = get_cpvtable_row(t, overwrite=overwrite)
        rows.append(r)

    rdf = pd.concat(rows)

    df = _df.merge(rdf, how='inner', on='ticid')
    df = df.reset_index(drop=True)

    msg = 'every unique ticid should have a row entry'
    assert len(_df) == len(df), msg

    csvpath = join(indir, "20230613_LGB_RJ_CPV_TABLE_supplemental.csv")
    df.to_csv(csvpath, index=False, sep="|")
    print(f"Wrote {csvpath}")

    # re-apply selection function
    sel = (
        (df.tic8_Tmag < 16)
        &
        (df.bp_rp > 1.5)
        &
        (df.M_G > 4)
        &
        (df.dist_pc < 150)
    )

    min_goodsector = (
        df.goodsectors.str.split(",").apply(
            lambda x: np.array(x).astype(int)
            if not np.all(pd.isnull(x)) else np.array([np.nan])
        ).apply(lambda x: np.nanmin(x))
    )
    min_maybesector = (
        df.maybesectors.str.split(",").apply(
            lambda x: np.array(x).astype(int)
            if not np.all(pd.isnull(x)) else np.array([np.nan])
        ).apply(lambda x: np.nanmin(x))
    )
    sel_max55 = (
        (min_goodsector <= 55) | (min_maybesector <= 55)
    )
    selcols = "ticid goodsectors maybesectors".split()
    print(f"\n...\n")
    print(f"WRN! including\n{df[~sel_max55][selcols]}\n, bc alpha per")
    print(f"\n...\n")

    N_0 = len(df)
    sdf = df[sel]

    sdf['qual'] = ~pd.isnull(sdf.goodsectors)
    sdf['quality'] = sdf['qual'].astype(int)

    # if "debunked" gets quality -1
    debunked_sel = sdf.ticid.isin(debunked_df.ticid)
    sdf.loc[debunked_sel, 'quality'] = -1


    csvpath = join(indir, "20230613_LGB_RJ_CPV_TABLE_supplemental_selfnapplied.csv")
    sdf.to_csv(csvpath, index=False, sep="|")
    print(f"Wrote {csvpath}")

    N_allcands = len(sdf)
    N_cqvs_nodebunked = len(sdf[(sdf.quality == 1) | (sdf.quality == 0)])

    N_2 = len(sdf[sdf.quality == 1])
    N_3 = len(sdf[sdf.quality == 0])
    N_debunked = len(sdf[sdf.quality == -1])

    # "good" CPVs

    sgdf = sdf[sdf.quality == 1]
    # "maybe" CPVs
    smdf = sdf[sdf.quality == 0]
    # debunked CPVs
    sddf = sdf[sdf.quality == -1]

    def num2word(num):
        units = [
            "zero", "one", "two", "three", "four", "five", "six",
            "seven", "eight", "nine", "ten"
        ]
        assert num <= 10
        return units[num]

    N_4 = num2word(len(sgdf[sgdf.banyan_assoc == 'FIELD']) + 1)
    N_5 = len(smdf[smdf.banyan_assoc == 'FIELD'])

    N_6 = len(sgdf[sgdf.banyan_assoc != 'FIELD']) - 1
    N_7 = len(smdf[smdf.banyan_assoc != 'FIELD'])

    N_8 = len(sgdf[sgdf.ruwe > 2])
    N_9 = len(smdf[smdf.ruwe > 2])

    _mdf = pd.concat((sgdf, smdf))
    dipcount_ticid = set(np.array(dipcount_df[dipcount_df.ticid.isin(_mdf.ticid)].ticid))
    fourier_ticid = set(np.array(fourier_df[fourier_df.ticid.isin(_mdf.ticid)].ticid))
    from gyrojo import venn
    plt.close("all")
    labels = venn.get_labels([dipcount_ticid, fourier_ticid],
                             fill=['number','logic'])
    fig,ax = venn.venn2(labels, names=['dipcount','fourier'])
    _outdir = join(RESULTSDIR, 'venn')
    if not os.path.exists(_outdir): os.mkdir(_outdir)
    _outpath = join(_outdir, 'dipcount_fourier_venn.png')
    fig.savefig(_outpath, bbox_inches='tight')
    print(f'wrote {_outpath}')
    plt.close("all")

    indip = _mdf.ticid.isin(dipcount_df.ticid)
    infourier = _mdf.ticid.isin(fourier_df.ticid)

    N_both = np.sum(indip & infourier)
    N_indip_notfourier = np.sum((indip) & (~infourier))
    N_infourier_notdip = np.sum((~indip) & (infourier))

    ticid_both = _mdf[indip & infourier].ticid
    print(f'ticid_both\n{ticid_both}')
    ticid_indip_notfourier = _mdf[(indip) & (~infourier)].ticid
    print(f'ticid_indip_notfourier\n{ticid_indip_notfourier}')
    ticid_infourier_notdip = _mdf[(~indip) & (infourier)].ticid
    print(f'ticid_infourier_notdip\n{ticid_infourier_notdip}')

    txt = (
        f"\n...\n"
        f"N={N_0} unique TICID's labelled 'good' or 'maybe' CPVs before selection fn"
        f"\n...\n"
        f"N={N_cqvs_nodebunked} unique TICID's labelled 'good' or 'maybe' CPVs after selection fn"
        f"\n...\n"
        f"N={N_2} unique TICID's labelled 'good' CPVs after selection fn"
        f"\n...\n"
        f"N={N_3} unique TICID's labelled 'maybe' CPVs after selection fn"
        f"\n...\n"
        f"N={N_4} good CPVs with 'field' designation from BANYAN"
        f"\n...\n"
        f"N={N_5} maybe CPVs with 'field' designation from BANYAN"
        f"\n...\n"
        f"N={N_6} good CPVs without 'field' designation from BANYAN"
        f"\n...\n"
        f"N={N_7} maybe CPVs without 'field' designation from BANYAN"
    )
    print(txt)

    latex_txt = (
        r"\newcommand{\nallcands}{"+str(N_allcands)+"}\n"
        r"\newcommand{\ncqvsnodebunked}{"+str(N_cqvs_nodebunked)+"}\n"
        r"\newcommand{\ngoods}{"+str(N_2)+"}\n"
        r"\newcommand{\nmaybes}{"+str(N_3)+"}\n"
        r"\newcommand{\ndebunked}{"+str(N_debunked)+"}\n"
        r"\newcommand{\ngoodsfieldbanyan}{"+str(N_4)+"}\n"
        r"\newcommand{\nmaybesfieldbanyan}{"+str(N_5)+"}\n"
        r"\newcommand{\ngoodsnotfieldbanyan}{"+str(N_6)+"}\n"
        r"\newcommand{\nmaybesnotfieldbanyan}{"+str(N_7)+"}\n"
        r"\newcommand{\nnotfieldbanyan}{"+str(N_6+N_7)+"}\n"
        r"\newcommand{\ngoodhighruwe}{"+str(N_8)+"}\n"
        r"\newcommand{\nmaybehighruwe}{"+str(N_9)+"}\n"
        r"\newcommand{\nbothdipfourier}{"+str(N_both)+"}\n"
        r"\newcommand{\nyesdipnofourier}{"+str(N_indip_notfourier)+"}\n"
        r"\newcommand{\nyesfouriernodip}{"+str(N_infourier_notdip)+"}\n"
    )

    outpath = join(PAPERDIR, 'vals_cqv_table_statistics.tex')
    with open(outpath, 'w') as f:
        f.writelines(latex_txt)
    print(f"Wrote {outpath}")

    # write portion that will go to latex
    sdf['period'] = sdf['tlc_mean_period']*24
    sdf['age'] = sdf['banyan_adopted_age']
    sdf['assoc'] = sdf['banyan_adopted_assoc']
    sdf['p_assoc'] = sdf['banyan_adopted_prob']
    shortcols = (
        "ticid "
        #"dr2_source_id ra dec "
        "tic8_Tmag dist_pc "
        "bp_rp ruwe period assoc "
        #"banyan_prob"
        "age teff_sedfit "
        "rstar_sedfit mstar_parsec Rcr_over_Rstar  "
        "quality N_sectors".split()
    )
    longcols = (
        "ticid "
        "dr2_source_id ra dec "
        "tic8_Tmag dist_pc "
        "bp_rp ruwe period assoc "
        #"banyan_prob "
        "age "
        "teff_sedfit teff_sedfit_perr teff_sedfit_merr "
        "rstar_sedfit rstar_sedfit_perr rstar_sedfit_merr "
        "mstar_parsec Rcr_over_Rstar  "
        "quality N_sectors".split()
    )
    pcols = shortcols + ['dist_metric_parsec']

    wdf = sdf[shortcols]
    wldf = sdf[longcols]
    pdf = sdf[pcols]

    wdf = wdf.sort_values(by=['quality','tic8_Tmag'], ascending=[False,True])
    wldf = wldf.sort_values(by=['quality','tic8_Tmag'], ascending=[False,True])
    pdf = pdf.sort_values(by=['quality','tic8_Tmag'], ascending=[False,True])

    rounddict = {
        'ra': 5,
        'dec': 5,
        'tic8_Tmag': 2,
        'ruwe': 2,
        'bp_rp': 3,
        'dist_pc': 1,
        'period': 2,
        'rstar_sedfit': 2,
        'mstar_parsec': 2,
        'Rcr_over_Rstar': 2
    }
    formatters = {}
    for k,v in rounddict.items():
        formatters[k] = lambda x: np.round(x, v)
    formatters['age'] = lambda x: int(x) if ~pd.isnull(x) else 'NaN'

    wdf = wdf.round(rounddict)
    wldf = wldf.round(rounddict)
    pdf = pdf.round(rounddict)

    csvpath = join(indir, "20230613_LGB_RJ_CPV_TABLE_selfnapplied_rounded.csv")
    wdf.to_csv(csvpath, index=False, sep="|")
    print(f"Wrote {csvpath}")

    texpath = join(PAPERDIR, "20230613_LGB_RJ_CPV_TABLE_selfnapplied_rounded_data.tex")
    wdf.to_latex(texpath, index=False, longtable=False, formatters=formatters)
    with open(texpath, 'r') as f:
        texlines = f.readlines()
    texlines = texlines[4:-2] # trim header and footer
    with open(texpath, 'w') as f:
        f.writelines(texlines)
    print(f"Wrote {texpath}")

    csvpath = join(indir, "20230613_LGB_RJ_CPV_TABLE_selfnapplied_rounded_longMRT.csv")
    wldf.to_csv(csvpath, index=False, sep="|")
    print(f"Wrote {csvpath}")

    pd.options.display.max_rows = 5000
    print(pdf)


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

    dr2_source_id = str(gdr2_df.dr2_source_id.iloc[0])

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
    agefloatvals = ",".join(list(sdf.agefloatval.astype(str)))
    age_refs = ",".join(list(sdf.age_reference.astype(str)))

    if len(assocs.split(',')) == 1 and assocs.split(',')[0] == 'FIELD':
        adopted_assoc = 'FIELD'
        adopted_age_float = np.nan
        adopted_prob = 1
        adopted_age_ref = ''

    else:

        # first listed association is field -> take second assoc
        if assocs.split(',')[0] == 'FIELD':
            adopted_assoc = assocs.split(",")[1]
            adopted_prob = probs.split(",")[1]
            adopted_age_float = agefloatvals.split(",")[1]
            adopted_age_ref = age_refs.split(',')[1]

        if assocs.split(',')[0] != 'FIELD':

            if len(assocs.split(',')) > 1:
                adopted_assoc = assocs.split(",")[0]
                adopted_prob = probs.split(",")[0]
                adopted_age_float = agefloatvals.split(",")[0]
                adopted_age_ref = age_refs.split(',')[0]

            else:
                adopted_assoc = assocs
                adopted_prob = probs
                adopted_age_float = agefloatvals
                adopted_age_ref = age_refs

    # alpha per case
    if dr2_source_id == '438223592047536640':
        assocs = 'APER'
        adopted_assoc = 'APER'
        banyan_prob = 0
        adopted_prob = 1
        ages = '$86 \pm 14$'
        adopted_age_float = 86
        age_refs = 'BoyleBouma2023'
        adopted_age_ref = 'BoyleBouma2023'

    bdf = pd.DataFrame({
        'banyan_assoc': assocs,
        'banyan_adopted_assoc': adopted_assoc,
        'banyan_prob': probs,
        'banyan_adopted_prob': adopted_prob,
        'banyan_age': ages,
        'banyan_adopted_age': float(adopted_age_float),
        'banyan_log_adopted_age': np.log10(float(adopted_age_float)*1e6),
        'banyan_age_refs': age_refs,
        'banyan_adopted_age_ref': adopted_age_ref,
    }, index=[0])

    return bdf


def get_tic8_row(t):

    cachepath = join(indir, f"TIC_{t}_mast_tic8_query.csv")

    if os.path.exists(cachepath):
        print(f'Found {cachepath}, returning.')
        return pd.read_csv(cachepath)

    ticstr = f"TIC {t}"
    MAX_ITER = 10
    ix = 0
    # MAST hosting anything is a recipe for failed queries, TBH.
    tic_data = None
    while ix < MAX_ITER and tic_data is None:
        try:
            tic_data = Catalogs.query_object(ticstr, catalog="TIC")
        except ResolverError as e:
            time.sleep(3)
            ix += 1
            print(f'TIC {t} failed initial MAST query with {e}, retrying {ix}/{MAX_ITER}...')
    if tic_data is None:
        raise ResolverError(f'TIC {t} failed to get MAST query.')

    t8_row = pd.DataFrame(tic_data.to_pandas().iloc[0]).T

    t8_row = t8_row.rename({c:f"tic8_{c}" for c in t8_row.columns},
                           axis='columns')

    t8_row.to_csv(cachepath, index=False)
    print(f'Cached {cachepath}')

    return t8_row


def identify_outliers(array, threshold=0.8):
    """
    Identifies outliers in an array of floats using a majority
    consensus style veto based on z-scores.

    Args:
        array (numpy.ndarray or list): The input array of floats.

        threshold (float, optional): The threshold for outlier
        identification. Defaults to 0.8.

    Returns:
        numpy.ndarray: A boolean mask indicating the outliers in the
        input array.
    """
    array = np.array(array)
    mean = np.mean(array)
    std = np.std(array)
    z_scores = (array - mean) / std

    consensus_threshold = int(threshold * len(array))
    outlier_mask = np.abs(z_scores) > 1

    while np.sum(outlier_mask) > consensus_threshold:
        array = array[~outlier_mask]
        mean = np.mean(array)
        std = np.std(array)
        z_scores = (array - mean) / std
        outlier_mask = np.abs(z_scores) > 1

    final_mask = np.zeros_like(z_scores, dtype=bool)
    final_mask[np.abs(z_scores) > 1] = True

    return final_mask



def get_tess_cpv_lc_properties(ticid):

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

        # NOTE: the properties dict has a bunch of information, including dip
        # depth.  this might be useful some day, but not today.
        #array = np.array
        #properties = eval(stc['properties'])

        nsplines_singlephase = stc['nsplines_singlephase']

        sectors.append(sector)
        periods.append(period)
        pgconds.append(periodogram_condition)
        n_peaks.append(_n_peaks)
        peak_phaseunits.append(_peaks_phaseunits)
        #peak_props.append(properties)
        nspls.append(nsplines_singlephase)

    periods = np.array(periods).astype(float)

    # ticids which upon manual inspection of the tlc_periods column will work
    # with identify_outliers this is the better method in that it's still
    # averaging over multiple sectors (whereas the choices made below, except
    # for TIC 402980664, are adopted from ONE sector)
    BAD_TICIDS_AUTO = [
        '177309964', '193136669', '2234692', '264767454', '300651846', '440725886', '57830249'
    ]
    #ticids which need the period set manually
    BAD_TICIDS_MANUAL = {
        '234295610': 0.76202,
        '272248916': 0.37084,
        '289840926': 0.199995,
        '335598085': 0.66034,
        '402980664': 18.5611/24,
    }

    use_periods = periods * 1.

    if str(ticid) in BAD_TICIDS_MANUAL:
        use_periods = np.array([ BAD_TICIDS_MANUAL[str(ticid)] ])
        msg = f'WRN! TIC {ticid} got {periods}, using {use_periods}'

    elif str(ticid) in BAD_TICIDS_AUTO:
        outliers = identify_outliers(periods)
        use_periods = periods[~outliers]
        msg = f'WRN! TIC {ticid} got {periods}, dropping {periods[outliers]}'
        print(42*'!')
        print(msg)
        print(42*'!')

    tlc_df = pd.DataFrame({
        'tlc_mean_period': np.nanmean(use_periods).round(6),
        'tlc_mean_npeaks': np.nanmean(n_peaks).round(1),
        'tlc_sectors': ",".join(np.array(sectors).astype(str)),
        'tlc_npeaks': ",".join(np.array(n_peaks).astype(str)),
        'tlc_periods': ",".join(np.array(periods).round(5).astype(str)),
        'tlc_useperiods': ",".join(np.array(use_periods).round(5).astype(str)),
        'tlc_peakloc_phaseunits': repr(peak_phaseunits).replace("],", "];")[1:-1]
    }, index=[0])

    return tlc_df


def get_sedfit_results(ticid):
    """
    teff_sedfit rstar_sedfit temp_suppression radius_inflation
    teff_noactivity rstar_noactivity
    """

    sed_dir = join(RESULTSDIR, "ariadne_sed_fitting", f"TIC_{ticid}")
    sed_path = join(sed_dir, "best_fit_average.dat")

    if not os.path.exists(sed_path):

        print(f"No SED fit result for {ticid} b/c ariadne has not run.")

        sed_df = pd.DataFrame({
            'teff_sedfit': 'PENDING',
            'teff_sedfit_perr': 'PENDING',
            'teff_sedfit_merr': 'PENDING',
            'rstar_sedfit': 'PENDING',
            'rstar_sedfit_perr':'PENDING',
            'rstar_sedfit_merr':'PENDING',
            'Av_sedfit': 'PENDING',
            'Av_sedfit_perr': 'PENDING',
            'Av_sedfit_merr': 'PENDING',
        }, index=[0])

        return sed_df

    else:

        a_df = pd.read_csv(sed_path, delim_whitespace=True)
        a_df.index = a_df['#Parameter']

    r = lambda x: int(np.round(x))
    r1 = lambda x: np.round(x, 3)

    sed_df = pd.DataFrame({
        'teff_sedfit': r(a_df.loc['teff', 'median']),
        'teff_sedfit_perr': r(a_df.loc['teff', 'upper']),
        'teff_sedfit_merr': r(a_df.loc['teff', 'lower']),
        'rstar_sedfit': r1(a_df.loc['rad', 'median']),
        'rstar_sedfit_perr': r1(a_df.loc['rad', 'upper']),
        'rstar_sedfit_merr': r1(a_df.loc['rad', 'lower']),
        'Av_sedfit': r1(a_df.loc['Av', 'median']),
        'Av_sedfit_perr': r1(a_df.loc['Av', 'upper']),
        'Av_sedfit_merr':  r1(a_df.loc['Av', 'lower']),
    }, index=[0])

    return sed_df


def get_isochrone_mass(sed_df, bdf):
    """
    Through naive nearest-neighbor interpolation against the PARSEC v1.2
    isochrones, guess the isochronal mass.

    The PARSEC grid was downloaded at a spacing of 1 to 200 myr, linearly
    spaced by 1 myr; all available masses.  Solar metallicity.
    """

    CONDITIONS = (
        (bdf['banyan_adopted_age'].iloc[0] == "???")
        or
        (pd.isnull(float(bdf['banyan_adopted_age'].iloc[0])))
        or
        (sed_df['rstar_sedfit'].iloc[0] == "PENDING")
    )
    if CONDITIONS:
        iso_df = pd.DataFrame({
            'mstar_parsec': np.nan,
            'logg_parsec': np.nan,
            'rstar_parsec': np.nan,
            'age_parsec': np.nan,
            'teff_parsec': np.nan,
            'dist_metric_parsec': np.nan,
        }, index=[0])
        return iso_df

    rstar = float(sed_df['rstar_sedfit'])
    rstar_err = np.mean([
        float(sed_df['rstar_sedfit_perr']),
        float(sed_df['rstar_sedfit_merr'])
    ])

    teff =  float(sed_df['teff_sedfit'])
    teff_err = np.mean([
        float(sed_df['teff_sedfit_perr']),
        float(sed_df['teff_sedfit_merr'])
    ])

    age = float(bdf['banyan_adopted_age'].iloc[0])
    if age < 10:
        # set age floor of 10 myr
        age = 10
    # require age uncertainty of 10% or +/-5 myr, whichever is bigger
    age_err = np.max([age*0.1, 5])

    df = get_PARSEC()

    dist = np.sqrt(
        np.array( ( (df['Rstar'] - rstar) / rstar_err )**2 )
        +
        np.array( ( (df['Teff'] - teff) / teff_err )**2 )
        +
        np.array( ( (df['age'] - age) / age_err )**2 )
    )

    df['dist'] = dist

    mstar_parsec = float(df.sort_values(by='dist').head(n=1)['Mass'])
    logg_parsec = float(df.sort_values(by='dist').head(n=1)['logg'])
    rstar_parsec = float(df.sort_values(by='dist').head(n=1)['Rstar'])
    age_parsec = float(df.sort_values(by='dist').head(n=1)['age'])
    teff_parsec = float(df.sort_values(by='dist').head(n=1)['Teff'])
    dist_metric = float(df.sort_values(by='dist').head(n=1)['dist'])

    iso_df = pd.DataFrame({
        'mstar_parsec': mstar_parsec,
        'logg_parsec': logg_parsec,
        'rstar_parsec': rstar_parsec,
        'age_parsec': age_parsec,
        'teff_parsec': teff_parsec,
        'dist_metric_parsec': dist_metric,
    }, index=[0])

    return iso_df


def get_cpvtable_row(ticid, overwrite=0):
    """
    ticid: str, e.g., "402980664"
    """

    cachecsv = join(indir, f"TIC{ticid}_cpvtable_row.csv")
    if os.path.exists(cachecsv) and not overwrite:
        return pd.read_csv(cachecsv, sep="|")

    gdr2_df = get_gaia_rows(ticid, allcols=1)

    tdf = assess_tess_holdings(ticid, outdir=indir)

    ftdf = flatten_tdf(tdf, ticid)

    bdf = get_banyan_result(gdr2_df)

    t8_df = get_tic8_row(ticid)

    tlc_df = get_tess_cpv_lc_properties(ticid)

    sed_df = get_sedfit_results(ticid)

    iso_df = get_isochrone_mass(sed_df, bdf)

    pd.options.display.max_rows = 5000

    row = pd.concat((ftdf, gdr2_df, bdf, t8_df, tlc_df, sed_df, iso_df), axis='columns')

    from astropy import constants as const, units as u
    omega = 2*np.pi / (float(row['tlc_mean_period'])*u.day)
    row['Rcr'] = (
        (const.G * float(row['mstar_parsec'])*u.Msun / omega**2)**(1/3)
    ).to(u.Rsun).value
    row['Rcr_over_Rstar'] = row['Rcr'] / row['rstar_sedfit']

    # TODO
    # TODO
    # TODO : need to add these
    # TODO
    # TODO
    bonus_cols = (
        "t_cross s_max manual_notes R10k_spectra_available"
    ).split()
    for c in bonus_cols:
        row[c] = 'PENDING'

    row.to_csv(cachecsv, index=False, sep="|")
    print(f"Wrote {cachecsv}")

    return row


if __name__ == "__main__":
    main(overwrite=0)
