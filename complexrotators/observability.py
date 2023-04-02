"""
Given tic ids, get some basic features of the stars.  (Nominally, the features
that we would care about in order to decide on ground-based follow-up
priorities)

get_gaia_rows
get_tess_stats
    given_ticid_get_variability_params
    given_ticid_get_period
check_tesspoint
check_astroplan_months_observable
get_bestmonth_hoursobservable
merge_to_observability_table
"""
import pandas as pd, numpy as np, matplotlib.pyplot as plt
import os, pickle
from glob import glob
from functools import partial
from os.path import join

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord


from astrobase.lcmodels.eclipses import invgauss_eclipses_curvefit_func
from astrobase.services.identifiers import tic_to_gaiadr2

from scipy.optimize import curve_fit

from cdips.utils.gaiaqueries import (
    given_source_ids_get_gaia_data,
    parallax_to_distance_highsn
)

from complexrotators.paths import RESULTSDIR, TARGETSDIR, TABLEDIR
from complexrotators import pipeline_utils as pu

def get_gaia_rows(ticid):
    """
    most important for observability -- gaia:
      * brightness
      * distance
      * color
    """

    dr2_source_id = tic_to_gaiadr2(ticid)

    runid = f"dr2_{dr2_source_id}"

    dr2_source_ids = np.array([np.int64(dr2_source_id)])
    gaia_r = given_source_ids_get_gaia_data(
        dr2_source_ids, runid, n_max=5, overwrite=False,
        enforce_all_sourceids_viable=True, savstr='', which_columns='*',
        table_name='gaia_source', gaia_datarelease='gaiadr2', getdr2ruwe=False
    )
    gdf_ruwe = given_source_ids_get_gaia_data(
        dr2_source_ids, runid+"_ruwe", n_max=5, overwrite=False,
        enforce_all_sourceids_viable=True, savstr='', which_columns='*',
        table_name='gaia_source', gaia_datarelease='gaiadr2', getdr2ruwe=True
    )

    SELCOLS = (
        'source_id,ra,dec,parallax,parallax_error,pmra,pmdec,'+
        'phot_g_mean_mag,phot_rp_mean_mag,phot_bp_mean_mag,bp_rp,'+
        'radial_velocity'
    ).split(',')

    outdf = gaia_r[SELCOLS]

    outdf['ruwe'] = float(gdf_ruwe['ruwe'])

    d_pc, upper_unc, lower_unc  = parallax_to_distance_highsn(
        float(outdf['parallax']),
        e_parallax_mas=float(outdf['parallax_error']),
        gaia_datarelease='gaia_dr2'
    )

    outdf['dist_pc'] = d_pc
    outdf['dist_pc_upper_unc'] = upper_unc
    outdf['dist_pc_lower_unc'] = lower_unc

    return outdf


def given_ticid_get_period(ticid):

    manualperiodfile = join(TARGETSDIR, 'ticids_manual_periods.csv')
    mpdf = pd.read_csv(manualperiodfile)

    # first, check for manually inserted periods
    row = mpdf.loc[mpdf.ticid.astype(str) == ticid]
    if len(row) > 0:
        period = float(row.manual_period)
        return period, 0.01

    # next, check the plot_22b_phase cache.  NOTE: by default, this is where
    # the meat is found.  otherwise, this getter would need to be reformatted
    # to pull the lightcurve and do that analysis from the light curves
    # themselves.

    check_pkl_cache = 0

    if check_pkl_cache:
        cachedir = join(RESULTSDIR, '*_phase', f'TIC_{ticid}')
        # all cadences, over all sectors
        pklpaths = glob(join(cachedir, f'*{ticid}*cr_periodsearch.pkl'))
        assert len(pklpaths) > 0

        pdicts = []
        for pklpath in pklpaths:
            with open(pklpath, 'rb') as f:
                pdict = pickle.load(f)
            pdicts.append(pdict)

        periods = [pdict['period'] for pdict in pdicts]

    check_log_cache = 1

    if check_log_cache:
        cachedir = join(RESULTSDIR, "cpvvetter", "lt_95pc_good_cpv_info")
        logpaths = glob(join(cachedir, f'*{ticid}*.log'))
        assert len(logpaths) > 0

        periods = []
        for logpath in logpaths:
            st = pu.load_status(logpath)
            if 'cpv_periodsearch_results' in st:
                this_period = float(st['cpv_periodsearch_results']['period'])
                periods.append(this_period)

    period = np.nanmedian(periods)
    period_stdev = np.nanstd(periods)

    return period, period_stdev


def given_ticid_get_variability_params(ticid, period_guess=None):

    # check the plot_22b_phase cache.  NOTE: by default, this is where the
    # content is found.  otherwise, this getter would need to be reformatted to
    # pull the lightcurve and do that analysis from the light curves
    # themselves.

    check_pkl_cache = 0

    if check_pkl_cache:

        cachedir = join(RESULTSDIR, '*_phase', f'TIC_{ticid}')
        # all cadences, over all sectors
        pklpaths = glob(join(cachedir, f'*{ticid}*dtr_lightcurve.pkl'))
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

    check_log_cache = 1

    if check_log_cache:
        cachedir = join(RESULTSDIR, "cpvvetter", "lt_95pc_good_cpv_info")
        logpaths = glob(join(cachedir, f'*{ticid}*.log'))
        assert len(logpaths) > 0, 'check cachedir...'

        a95_5s = []
        for logpath in logpaths:
            st = pu.load_status(logpath)
            if 'count_phased_local_minima_results' in st:
                this_a = float(st['count_phased_local_minima_results']['a_95_5'])
                a95_5s.append(this_a)
        a_5_95 = np.nanmedian(a95_5s)

    return a_5_95


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
    a_5_95 = given_ticid_get_variability_params(
        ticid, period_guess=period
    )

    outdf = pd.DataFrame({
        'period':period,
        'period_stdev':period_stdev,
        'a_5_95':a_5_95,
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


def check_astroplan_months_observable(
    starid, ra, dec,
    site = 'keck',
    min_altitude = 30*u.deg,
    twilight_limit = 'astronomical',
    minokmoonsep=30*u.deg
):
    """
    Given a star ID, ra/dec, and a site location, figure out the best
    months for observing a star.

    Args:
        starid: string, e.g. "TIC_12345678"
        ra/dec: float degrees
    """

    from astroplan import (FixedTarget, Observer, is_observable,
                           months_observable,
                           AtNightConstraint, AltitudeConstraint,
                           LocalTimeConstraint, MoonSeparationConstraint,
                           AirmassConstraint, moon)

    target_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    target = FixedTarget(coord=target_coord, name=starid)
    observer = Observer.at_site(site)

    if twilight_limit == 'astronomical':
        twilight_constraint = AtNightConstraint.twilight_astronomical()
    elif twilight_limit == 'nautical':
        twilight_constraint = AtNightConstraint.twilight_nautical()
    else:
        raise NotImplementedError('civil twilight is janky.')

    constraints = [twilight_constraint,
                   AltitudeConstraint(min=min_altitude),
                   MoonSeparationConstraint(min=minokmoonsep)]

    bestmonths = months_observable(constraints, observer, [target])

    # computed observability on "bestmonths" grid of 0.5 hr
    print(f'for {starid}, got best-months on 0.5 hour grid:')
    print(bestmonths)
    print('where 1 = Jan, 2 = Feb, etc.')

    joiner = lambda x: ','.join(np.array(x).astype(str))

    if bestmonths == [set()]:
        bestmonths = str(-1)
    else:
        bestmonths = joiner(list(bestmonths[0]))

    outdf = pd.DataFrame({
        'bestmonths': bestmonths,
    }, index=[0])

    return outdf


def get_bestmonth_hoursobservable(
    starid, ra, dec,
    site = 'keck',
    min_altitude = 35*u.deg,
    twilight_limit = 'astronomical',
    minokmoonsep=30*u.deg,
    semester='22B'
):
    """
    Given a star ID, ra/dec, an observing run start/stop time, and a site
    location, figure out how many HOURS the star will observable at
    >40 deg altitude per night.

    Args:
        starid: string, e.g. "TIC_12345678"
        ra/dec: float degrees
    """

    from astroplan import (FixedTarget, Observer, is_observable,
                           months_observable,
                           AtNightConstraint, AltitudeConstraint,
                           LocalTimeConstraint, MoonSeparationConstraint,
                           AirmassConstraint, moon,
                           observability_table)

    target_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    target = FixedTarget(coord=target_coord, name=starid)
    observer = Observer.at_site(site)

    if twilight_limit == 'astronomical':
        twilight_constraint = AtNightConstraint.twilight_astronomical()
    elif twilight_limit == 'nautical':
        twilight_constraint = AtNightConstraint.twilight_nautical()
    else:
        raise NotImplementedError('civil twilight is janky.')

    constraints = [twilight_constraint,
                   AltitudeConstraint(min=min_altitude),
                   MoonSeparationConstraint(min=minokmoonsep)]

    if semester.endswith("B"):
        # B semester example:
        #min_time = Time('2022-07-31 23:59:00'),
        #max_time = Time('2023-01-31 23:59:00')
        starttime = Time('2022-08-15 23:59:00')
        times = [starttime + ix*30*u.day for ix in range(6)]
    elif semester.endswith("A"):
        pass

    obs_tables = {}
    hrs_visible_per_night = {}
    for time in times:
        obs_table = observability_table(
            constraints,
            observer,
            [target],
            time_range=time,
            time_grid_resolution=0.25*u.hour
        )
        obs_tables[time] = obs_table
        # hours observable per night at this time
        month_key = time.to_datetime().month
        hrs_visible_per_night[month_key] = np.round(float(
            obs_table['fraction of time observable']
        )*24,1)

    outdf = pd.DataFrame(hrs_visible_per_night, index=[0])

    return outdf




def merge_to_observability_table(
    csvpath, r_gaia, r_tess, r_tic8, r_tesspoint, r_astroplan
    ):

    outdf = pd.concat(
        (r_gaia, r_tess, r_tic8, r_tesspoint, r_astroplan), axis=1
    )

    outdf.to_csv(csvpath, index=False)
    print(f'Wrote {csvpath}')
