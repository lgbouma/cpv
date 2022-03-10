"""
Given tic ids, get some basic features of the stars.  (Nominally, the features
that we would care about in order to decide on ground-based follow-up
priorities)

get_gaia_rows
get_tess_stats
    given_ticid_get_variability_params
    given_ticid_get_period
check_tesspoint
check_astroplan_observability
merge_to_obsevability_table
"""
import pandas as pd, numpy as np, matplotlib.pyplot as plt
import os, pickle
from glob import glob
from functools import partial

from astrobase.lcmodels.eclipses import invgauss_eclipses_curvefit_func
from astrobase.services.identifiers import tic_to_gaiadr2

from scipy.optimize import curve_fit

from cdips.utils.gaiaqueries import (
    given_source_ids_get_gaia_data,
    parallax_to_distance_highsn
)

from complexrotators.paths import RESULTSDIR, TARGETSDIR, TABLEDIR

def get_gaia_rows(ticid):
    """
    most important for observability -- gaia:
      * brightness
      * distance
      * color
    """

    source_id = tic_to_gaiadr2(ticid)

    groupname = f'gaia_dr2_{source_id}'
    source_ids = np.array([np.int64(source_id)])

    gaia_r = given_source_ids_get_gaia_data(
        source_ids, groupname, n_max=10000, overwrite=False,
        enforce_all_sourceids_viable=True, savstr='', whichcolumns='*',
        gaia_datarelease='gaiaedr3', getdr2ruwe=False
    )

    SELCOLS = (
        'source_id,ra,dec,parallax,parallax_error,pmra,pmdec,'+
        'phot_g_mean_mag,phot_rp_mean_mag,phot_bp_mean_mag,bp_rp,'+
        'dr2_radial_velocity,ruwe'
    ).split(',')

    outdf = gaia_r[SELCOLS]

    d_pc, upper_unc, lower_unc  = parallax_to_distance_highsn(
        float(outdf['parallax']),
        e_parallax_mas=float(outdf['parallax_error']),
        gaia_datarelease='gaia_edr3'
    )

    outdf['dist_pc'] = d_pc
    outdf['dist_pc_upper_unc'] = upper_unc
    outdf['dist_pc_lower_unc'] = lower_unc

    return outdf


def given_ticid_get_period(ticid):

    manualperiodfile = os.path.join(TARGETSDIR, 'ticids_manual_periods.csv')
    mpdf = pd.read_csv(manualperiodfile)

    # first, check for manually inserted periods
    row = mpdf.loc[mpdf.ticid.astype(str) == ticid]
    if len(row) > 0:
        period = float(row.manual_period)
        return period

    # next, check the plot_22b_phase cache.  NOTE: by default, this is where
    # the meat is found.  otherwise, this getter would need to be reformatted
    # to pull the lightcurve and do that analysis from the light curves
    # themselves.

    cachedir = os.path.join(RESULTSDIR, '22b_phase', f'TIC_{ticid}')
    # all cadences, over all sectors
    pklpaths = glob(os.path.join(cachedir, f'*{ticid}*cr_periodsearch.pkl'))
    assert len(pklpaths) > 0

    pdicts = []
    for pklpath in pklpaths:
        with open(pklpath, 'rb') as f:
            pdict = pickle.load(f)
        pdicts.append(pdict)

    periods = [pdict['period'] for pdict in pdicts]

    period = np.nanmedian(periods)
    period_stdev = np.nanstd(periods)

    return period, period_stdev


def given_ticid_get_variability_params(ticid, period_guess=None):

    # check the plot_22b_phase cache.  NOTE: by default, this is where the
    # content is found.  otherwise, this getter would need to be reformatted to
    # pull the lightcurve and do that analysis from the light curves
    # themselves.

    cachedir = os.path.join(RESULTSDIR, '22b_phase', f'TIC_{ticid}')
    # all cadences, over all sectors
    pklpaths = glob(os.path.join(cachedir, f'*{ticid}*dtr_lightcurve.pkl'))
    assert len(pklpaths) > 0

    lcdicts = []
    for pklpath in pklpaths:
        with open(pklpath, 'rb') as f:
            lcdict = pickle.load(f)
        lcdicts.append(lcdict)

    get_5_95 = lambda x: np.nanpercentile(x, 95) - np.nanpercentile(x, 5)
    get_10_90 = lambda x: np.nanpercentile(x, 90) - np.nanpercentile(x, 10)

    a_5_95 = np.nanmedian([get_5_95(lcdict['y_flat']) for lcdict in lcdicts])
    a_10_90 = np.nanmedian([get_10_90(lcdict['y_flat']) for lcdict in lcdicts])

    fit_param_list = []

    #TODO might want to implement like a double fit gaussian or something
    #for lcdict in lcdicts:
    #    time = lcdict['x_trend']
    #    flux = lcdict['y_flat']

    #    initial_params = {'period': period_guess,
    #                      'epoch': 1618,
    #                      'pdepth': 0.05,
    #                      'pduration': 0.05,
    #                      'psdepthratio': 2,
    #                      'secondaryphase': 0.5}

    #    fit_params, fit_cov = curve_fit(
    #        invgauss_eclipses_curvefit_func, time, flux-np.nanmedian(flux)
    #    )
    #    fit_param_list.append(fit_params)

    #import IPython; IPython.embed()
    #assert 0

    return a_5_95, a_10_90


def get_tess_stats(ticid):
    """
    most important for observability -- tess:
      * period
      * dip depth
      * dip duration(?)
      * number of dips per cycle(?)
    """

    # get the period
    period, period_stdev = given_ticid_get_period(ticid)

    # get the variability amplitude diagnostics
    a_5_95, a_10_90 = given_ticid_get_variability_params(
        ticid, period_guess=period
    )

    outdf = pd.DataFrame({
        'period':period,
        'period_stdev':period_stdev,
        'a_5_95':a_5_95,
        'a_10_90':a_10_90
    }, index=[0])

    return outdf


def check_tesspoint(ra, dec, ticid):

    from tess_stars2px import tess_stars2px_function_entry

    (outID, outEclipLong, outEclipLat,
     sector, cam, ccd,
     colpix, rowpix, scinfo ) = (
         tess_stars2px_function_entry(ticid, ra, dec)
     )

    joiner = lambda x: ','.join(np.array(x).astype(str))
    pxjoiner = lambda x: ','.join(np.array(np.round(x,1)).astype(str))

    outdf = pd.DataFrame({
        'sector': joiner(sector),
        'cam': joiner(cam),
        'ccd': joiner(ccd),
        'colpix': pxjoiner(colpix),
        'rowpix': pxjoiner(rowpix)
    }, index=[0])

    return outdf


def check_astroplan_observability(ticid):

    #FIXME FIXME in the morning...

    #TODO
    pass

def merge_to_obsevability_table(ticid):
    #TODO
    pass
