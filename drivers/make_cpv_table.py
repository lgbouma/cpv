"""
Construct the extended and publication versions of the 2023 catalog table.

Uses log output created by running find_CPVs.py using
sample_id=="2023catalog_LGB_RJ_concat" first (see HOWTO.md)

Utilities:

| get_cpvtable_row(ticid)
| gdr2_df = get_gaia_dr2_rows(ticid)
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

from astropy import constants as const, units as u

from complexrotators.observability import (
    get_gaia_dr2_rows, assess_tess_holdings, get_gaia_dr3_rows
)
from complexrotators.getters import get_tic8_row

from complexrotators import pipeline_utils as pu

vetdir = join(RESULTSDIR, "cpvvetter")
indir = join(TABLEDIR, "2023_catalog_table")
assert os.path.exists(indir)

def main(overwrite=0):

    csvpath = join(indir, "20230613_LGB_RJ_uticid_quality_label.csv")
    _df = pd.read_csv(csvpath, sep="|")

    csvpath2 = join(DATADIR, "targetlists", "20230813_DEBUNKED_CQVs.txt")
    debunked_df = pd.read_csv(csvpath2)

    csvpath3 = join(DATADIR, "targetlists",
                    "20230411_goodandmaybe_CPV_ticids_d_lt_150pc.csv")
    dipcount_df = pd.read_csv(csvpath3)

    csvpath4 = join(DATADIR, "targetlists",
                    "20230501_RAHUL_FULL_LIST_NO_DUPLICATES.csv")
    fourier_df = pd.read_csv(csvpath4)

    csvpath5 = join(TABLEDIR, "multiple_period_20230613_LGB_RJ_CPV",
                    "20230814_multiperiod_MANUAL.csv")
    multperiod_df = pd.read_csv(csvpath5)
    selcols = 'ticid P2_hr P2class'
    multperiod_df = multperiod_df[selcols.split()]

    HACK = 1
    if HACK: # once things have already been made...
        __df = pd.read_csv("table1_short.csv", sep="|")
        _df = _df[_df.ticid.isin(__df.ticid)]

    rows = []
    for t in _df['ticid']:
        r = get_cpvtable_row(t, overwrite=overwrite)
        rows.append(r)
        assert r.shape == (1, 325)

    rdf = pd.concat(rows).reset_index(drop=True)

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

    sel &= sel_max55

    selcols = "ticid goodsectors maybesectors".split()
    print(f"\n...\n")
    print(f"WRN! dropping {df[~sel_max55][selcols]}\n")
    print(f"\n...\n")

    N_0 = len(df)
    sdf = df[sel]

    sdf['qual'] = ~pd.isnull(sdf.goodsectors)
    sdf['quality'] = sdf['qual'].astype(int)

    # if "debunked" gets quality -1
    debunked_sel = sdf.ticid.isin(debunked_df.ticid)
    sdf.loc[debunked_sel, 'quality'] = -1

    # sort
    sdf = sdf.sort_values(by=['quality','tic8_Tmag'], ascending=[False,True])

    # update dr3 distance and RUWE for TIC 368129164 manually
    sdf.loc[sdf.ticid == 368129164, 'dr3_dist_pc'] = (
        sdf.loc[sdf.ticid == 368129164, 'dist_pc']
    )
    sdf.loc[sdf.ticid == 368129164, 'dr3_ruwe'] = (
        sdf.loc[sdf.ticid == 368129164, 'ruwe']
    )

    # assert no obvious binaries
    assert sdf.dr3_non_single_star.sum() == 0

    # merge for secondary periods
    N_before = len(sdf)
    sdf = sdf.merge(multperiod_df, how='inner', on='ticid')
    N_after = len(sdf)
    sdf = sdf.reset_index(drop=True)

    msg = 'every ticid should have a P2 entry, even if nan'
    assert N_before == N_after, msg


    # construct the binarity indicator
    rvscatterflag = (sdf.dr3_radial_velocity_error > 20)
    ruweflag = (sdf.dr3_ruwe > 2)
    weakruweflag = (sdf.dr3_ruwe > 1.4)
    multipleperiodflag = (
        sdf.ticid.isin(multperiod_df[multperiod_df.P2_hr != '-'].ticid)
    )

    binaryflag_strs = []
    sdf['rvscatterflag'] = rvscatterflag.astype(int).astype(str)
    sdf['ruweflag'] = ruweflag.astype(int).astype(str)
    sdf['weakruweflag'] = weakruweflag.astype(int).astype(str)
    sdf['multipleperiodflag'] = multipleperiodflag.astype(int).astype(str)

    # rv scatter, RUWE, multiple period
    sdf['binarityflag'] = (
        sdf['rvscatterflag'] + sdf['ruweflag']  + sdf['multipleperiodflag']
    )

    ## manual string casting
    #strcols = 'dr2_source_id dr3_source_id'.split()
    #for strcol in strcols:
    #    sdf[strcol] = '"' + sdf[strcol].astype(str) + '"'

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

    N_4 = num2word(len(sgdf[sgdf.banyan_assoc == 'FIELD']))
    N_5 = len(smdf[smdf.banyan_assoc == 'FIELD'])

    N_6 = len(sgdf[sgdf.banyan_assoc != 'FIELD'])
    N_7 = len(smdf[smdf.banyan_assoc != 'FIELD'])

    N_8 = len(sgdf[sgdf.dr3_ruwe > 2])
    N_9 = len(smdf[smdf.dr3_ruwe > 2])

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

    ticid_both = np.array(_mdf[indip & infourier].ticid)
    print(f'ticid_both\n{ticid_both}')
    ticid_indip_notfourier = np.array(_mdf[(indip) & (~infourier)].ticid)
    print(f'ticid_indip_notfourier\n{ticid_indip_notfourier}')
    ticid_infourier_notdip = np.array(_mdf[(~indip) & (infourier)].ticid)
    print(f'ticid_infourier_notdip\n{ticid_infourier_notdip}')

    #
    # binarity counts
    #
    N_rv_scatter = len(_mdf[_mdf.rvscatterflag.astype(int) == 1])
    N_ruwe = len(_mdf[_mdf.ruweflag.astype(int) == 1])
    N_weakruwe = len(_mdf[_mdf.weakruweflag.astype(int) == 1])
    N_multperiod = len(_mdf[_mdf.multipleperiodflag.astype(int) == 1])

    N_goodruweandmultperiod = len(_mdf[
        (_mdf.multipleperiodflag.astype(int) == 1)
        &
        (_mdf.ruweflag.astype(int) == 1)
        &
        (_mdf.quality == 1)
    ])
    N_mayberuweandmultperiod = len(_mdf[
        (_mdf.multipleperiodflag.astype(int) == 1)
        &
        (_mdf.ruweflag.astype(int) == 1)
        &
        (_mdf.quality == 0)
    ])
    N_goodweakruwe = len(_mdf[
        (_mdf.weakruweflag.astype(int) == 1)
        &
        (_mdf.quality == 1)
    ])
    N_maybeweakruwe = len(_mdf[
        (_mdf.weakruweflag.astype(int) == 1)
        &
        (_mdf.quality == 0)
    ])
    N_goodweakruweandmultperiod = len(_mdf[
        (_mdf.multipleperiodflag.astype(int) == 1)
        &
        (_mdf.weakruweflag.astype(int) == 1)
        &
        (_mdf.quality == 1)
    ])
    N_maybeweakruweandmultperiod = len(_mdf[
        (_mdf.multipleperiodflag.astype(int) == 1)
        &
        (_mdf.weakruweflag.astype(int) == 1)
        &
        (_mdf.quality == 0)
    ])

    N_goodmultperiod = len(_mdf[
        (_mdf.multipleperiodflag.astype(int) == 1)
        &
        (_mdf.quality == 1)
    ])
    N_maybemultperiod = len(_mdf[
        (_mdf.multipleperiodflag.astype(int) == 1)
        &
        (_mdf.quality == 0)
    ])


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
        r"\newcommand{\ngoodweakruwe}{"+str(N_goodweakruwe)+"}\n"
        r"\newcommand{\nmaybeweakruwe}{"+str(N_maybeweakruwe)+"}\n"
        r"\newcommand{\nbothdipfourier}{"+str(N_both)+"}\n"
        r"\newcommand{\nyesdipnofourier}{"+str(N_indip_notfourier)+"}\n"
        r"\newcommand{\nyesfouriernodip}{"+str(N_infourier_notdip)+"}\n"
        r"\newcommand{\nrvscatterflag}{"+str(num2word(N_rv_scatter))+"}\n"
        r"\newcommand{\nruweflag}{"+str(N_ruwe)+"}\n"
        r"\newcommand{\nweakruweflag}{"+str(N_weakruwe)+"}\n"
        r"\newcommand{\nmultperiodflag}{"+str(N_multperiod)+"}\n"
        r"\newcommand{\ngoodmultperiodflag}{"+str(N_goodmultperiod)+"}\n"
        r"\newcommand{\nmaybemultperiodflag}{"+str(N_maybemultperiod)+"}\n"
        r"\newcommand{\ngoodruweandmultperiod}{"+str(N_goodruweandmultperiod)+"}\n"
        r"\newcommand{\nmayberuweandmultperiod}{"+str(N_mayberuweandmultperiod)+"}\n"
        r"\newcommand{\ngoodweakruweandmultperiod}{"+str(N_goodweakruweandmultperiod)+"}\n"
        r"\newcommand{\nmaybeweakruweandmultperiod}{"+str(N_maybeweakruweandmultperiod)+"}\n"
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
        "tic8_Tmag dr3_dist_pc "
        "dr3_bp_rp dr3_ruwe period assoc "
        #"banyan_prob"
        "age teff_sedfit "
        "rstar_sedfit mass_parsec Rcr_over_Rstar  "
        "P2_hr "
        "quality binarityflag N_sectors".split()
    )
    longcols = (
        "ticid "
        "dr2_source_id dr3_source_id angular_distance magnitude_difference dr3_ra dr3_dec "
        "tic8_Tmag dr3_dist_pc "
        "dr3_bp_rp dr3_ruwe period assoc "
        #"banyan_prob "
        "age "
        "teff_sedfit teff_sedfit_perr teff_sedfit_merr "
        "rstar_sedfit rstar_sedfit_perr rstar_sedfit_merr "
        "mass_parsec mass_parsec_perr mass_parsec_merr "
        "Rcr_over_Rstar "
        "P2_hr "
        "quality binarityflag N_sectors".split()
    )
    pcols = shortcols + ['dist_median_parsec']

    wdf = sdf[shortcols]
    wldf = sdf[longcols]
    pdf = sdf[pcols]

    wdf = wdf.sort_values(by=['quality','tic8_Tmag'], ascending=[False,True])
    wldf = wldf.sort_values(by=['quality','tic8_Tmag'], ascending=[False,True])
    pdf = pdf.sort_values(by=['quality','tic8_Tmag'], ascending=[False,True])

    rounddict = {
        'dr3_ra': 5,
        'dr3_dec': 5,
        'tic8_Tmag': 2,
        'dr3_ruwe': 2,
        'dr3_bp_rp': 3,
        'dist_pc': 1,
        'dr3_dist_pc': 1,
        'period': 2,
        'rstar_sedfit': 2,
        'rstar_sedfit_perr': 2,
        'rstar_sedfit_merr': 2,
        'mass_parsec': 2,
        'mass_parsec_perr': 2,
        'mass_parsec_merr': 2,
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

    selcols = (
        'ticid rstar_sedfit rstar_sedfit_perr rstar_sedfit_merr '
        'teff_sedfit teff_sedfit_perr teff_sedfit_merr '
        'mass_parsec mass_parsec_perr mass_parsec_merr '
    ).split()
    print(42*'-')
    print(wldf[selcols])

    teff_err = np.nanmean(np.vstack([
        np.array(wldf['teff_sedfit_perr']),
        np.array(wldf['teff_sedfit_merr'])
    ]), axis=0)
    teff_err_rel = np.median( 100 * teff_err / np.array(wldf['teff_sedfit']) )

    rstar_err = np.nanmean(np.vstack([
        np.array(wldf['rstar_sedfit_perr']),
        np.array(wldf['rstar_sedfit_merr'])
    ]), axis=0)
    rstar_err_rel = np.median( 100 * rstar_err / np.array(wldf['rstar_sedfit']))

    mass_err = np.nanmean(np.vstack([
        np.array(wldf['mass_parsec_perr']),
        np.array(wldf['mass_parsec_merr'])
    ]), axis=0)
    mass_err_rel = np.nanmedian( 100 * mass_err / np.array(wldf['mass_parsec']) )

    print(f'median teff err abs: {np.median(teff_err):.2f} K%')
    print(f'median teff err: {teff_err_rel:.2f} %')
    print(f'median rstar err: {rstar_err_rel:.2f} %')
    print(f'median mass err: {mass_err_rel:.2f} %')




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

    # TIC 397791443: BANYAN assigned LCC prob 0.0008, is in IC2602 from CantatGaudin2020
    if dr2_source_id == '5239135492911437568':
        assocs = 'IC2602'
        adopted_assoc = 'IC2602'
        banyan_prob = 0
        adopted_prob = 1
        ages = '$46 \pm 6$'
        adopted_age_float = 46
        age_refs = 'Dobbie2010'
        adopted_age_ref = 'Dobbie2010'

    if str(dr2_source_id) == '860453786736413568':
        assocs = 'FIELD'
        adopted_assoc = 'FIELD'
        banyan_prob = 0
        adopted_prob = 1
        ages = '$16^{+19}_{-6}$'
        adopted_age_float = 16
        age_refs = 'Bouma2025'
        adopted_age_ref = 'Bouma2025'

    # MANUAL CASES WITH UNLIKELY MEMBERSHIPS (<1% and couldn't manually verify)
    if dr2_source_id in ['3311867153305660288', '5295268619510973056']:
        adopted_assoc = adopted_assoc + "(?)"

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

    sample_id = "2023catalog_LGB_RJ_concat_BACKUP"
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


def get_sedfit_results(ticid, uniformprior=1):
    """
    teff_sedfit rstar_sedfit temp_suppression radius_inflation
    teff_noactivity rstar_noactivity
    """

    sed_dir = join(RESULTSDIR, "ariadne_sed_fitting", f"TIC_{ticid}")
    if uniformprior:
        sed_dir = join(RESULTSDIR, "ariadne_sed_fitting_UNIFORM", f"TIC_{ticid}")
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
            'mass_parsec': np.nan,
            'logg_parsec': np.nan,
            'rstar_parsec': np.nan,
            'age_parsec': np.nan,
            'teff_parsec': np.nan,
            'dist_median_parsec': np.nan,
        }, index=[0])
        params = 'Mass logg Rstar age Teff dist_median'.split()
        for p in params:
            iso_df[f'{p.lower()}_parsec_perr'] = np.nan
            iso_df[f'{p.lower()}_parsec_merr'] = np.nan
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

    from complexrotators.isochroneinterp import nn_PARSEC_interpolator

    paramdict = nn_PARSEC_interpolator(
        rstar, teff, age, rstar_err, teff_err, age_err
    )

    iso_df = pd.DataFrame({}, index=[0])
    params = 'Mass logg Rstar age Teff dist_median'.split()
    for p in params:
        iso_df[f'{p.lower()}_parsec'] = paramdict[p][0]
        iso_df[f'{p.lower()}_parsec_perr'] = paramdict[p][1]
        iso_df[f'{p.lower()}_parsec_merr'] = paramdict[p][2]

    return iso_df


def get_cpvtable_row(ticid, overwrite=0):
    """
    ticid: str, e.g., "402980664"
    """

    cachecsv = join(indir, f"TIC{ticid}_cpvtable_row.csv")
    if os.path.exists(cachecsv) and not overwrite:
        return pd.read_csv(cachecsv, sep="|")

    gdr2_df = get_gaia_dr2_rows(ticid, allcols=1)

    gdr3_df = get_gaia_dr3_rows(ticid)

    tdf = assess_tess_holdings(ticid, outdir=indir)

    ftdf = flatten_tdf(tdf, ticid)

    bdf = get_banyan_result(gdr2_df)

    t8_df = get_tic8_row(ticid, indir)

    tlc_df = get_tess_cpv_lc_properties(ticid)

    sed_df = get_sedfit_results(ticid, uniformprior=1)

    iso_df = get_isochrone_mass(sed_df, bdf)

    pd.options.display.max_rows = 5000

    row = pd.concat((ftdf, gdr2_df, gdr3_df, bdf, t8_df, tlc_df, sed_df, iso_df), axis='columns')

    omega = 2*np.pi / (float(row['tlc_mean_period'])*u.day)
    row['Rcr'] = (
        (const.G * float(row['mass_parsec'])*u.Msun / omega**2)**(1/3)
    ).to(u.Rsun).value

    mass_err = np.nanmean([float(row['mass_parsec_perr']),
                           float(row['mass_parsec_merr'])])
    row['Rcr_err'] = float(row['Rcr']) * (
        (1/3) * (mass_err / float(row['mass_parsec']) )
    )

    row['Rcr_over_Rstar'] = row['Rcr'] / row['rstar_sedfit']

    rstar_sedfit_err = np.nanmean([float(row['rstar_sedfit_perr']),
                                   float(row['rstar_sedfit_merr'])])
    row['Rcr_over_Rstar_err'] = float(row['Rcr_over_Rstar']) * np.sqrt(
        (row['Rcr_err'] / row['Rcr'])**2
        +
        (rstar_sedfit_err / row['rstar_sedfit'] )**2
    )

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

    if not row.shape == (1, 327):
        import IPython; IPython.embed()
        assert 0

    return row


if __name__ == "__main__":
    ticid = '141146667'
    r = get_cpvtable_row(ticid, overwrite=0)
    import IPython; IPython.embed()
    assert 0
    main(overwrite=0)
