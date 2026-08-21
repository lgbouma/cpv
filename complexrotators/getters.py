"""
Contents:

General:

    | get_tic8_row
    | get_2min_cadence_spoc_tess_lightcurve
    | get_20sec_cadence_spoc_tess_lightcurve

    | _get_local_lcpaths_given_ticid
    | _get_lcpaths_given_ticid

    | get_cqv_search_sample

tic4029 specialized:

    | get_tic4029_lc_and_mask
    | get_4029_manual_mask

Spectra:

    | get_specriver_data
    | get_bochanski2007_m_standard
"""
#############
## LOGGING ##
#############
import logging
from complexrotators import log_sub, log_fmt, log_date_fmt

LOCAL_DEBUG = 0
DEBUG = False
if DEBUG:
    level = logging.DEBUG
else:
    level = logging.INFO
LOGGER = logging.getLogger(__name__)
logging.basicConfig(
    level=level,
    style=log_sub,
    format=log_fmt,
    datefmt=log_date_fmt,
    force=True
)

LOGDEBUG = LOGGER.debug
LOGINFO = LOGGER.info
LOGWARNING = LOGGER.warning
LOGERROR = LOGGER.error
LOGEXCEPTION = LOGGER.exception

#############
## IMPORTS ##
#############

import numpy as np, pandas as pd
import lightkurve as lk
import os, csv
from typing import Tuple
from os.path import join
from glob import glob
import subprocess
from complexrotators.paths import (
    LKCACHEDIR, LOCALDIR, RESULTSDIR, DATADIR, TARSCACHEDIR
)

from astropy.io import fits
from astropy import units as u, constants as const
from astropy.time import Time
from scipy.ndimage import gaussian_filter1d

from astroquery.exceptions import ResolverError
from astroquery.mast import Catalogs

from cdips.utils.lcutils import astropy_utc_time_to_bjd_tdb

def get_cqv_search_sample():
    localdir = os.path.expanduser('~/local/SPOCLC')
    csvpath = join(localdir, 'gaia_X_spoc2min_merge.csv')
    df = pd.read_csv(csvpath)

    sel = (
        (df.TESSMAG < 16) & (df.bp_rp > 1.5) &
        (df.M_G > 4) & (df.parallax > 1e3*(1/150)) &
        (df.SECTOR <= 55)
    )

    sdf = df[sel]
    return sdf


def get_tic4029_lc_and_mask(model_id):
    """
    model_id, e.g., 'manual_20230617_mask_v0_nterms2'
    """

    from complexrotators.lcprocessing import prepare_cpv_light_curve

    ticid = '402980664'
    sample_id = '2023catalog_LGB_RJ_concat' # used for cacheing
    cachedir = join(LOCALDIR, "cpv_finding", sample_id)

    lcpaths = _get_lcpaths_fromlightkurve_given_ticid(ticid)

    # stack over all sectors
    times, fluxs, cadencenos = [], [], []
    for lcpath in np.sort(lcpaths):
        (time, flux, qual, x_obs, y_obs, y_flat, y_trend, x_trend, cadence_sec,
         sector, starid, cadenceno) = prepare_cpv_light_curve(
             lcpath, cachedir, returncadenceno=1
         )
        times.append(x_obs)
        fluxs.append(y_flat)
        cadencenos.append(cadenceno)
    times = np.hstack(np.array(times, dtype=object).flatten())
    fluxs = np.hstack(np.array(fluxs, dtype=object).flatten())
    cadencenos = np.hstack(np.array(cadencenos, dtype=object).flatten())

    if 'manual_20230617_mask_v0' in model_id:
        full_ood = get_4029_manual_mask(times, fluxs, cadencenos, get_full=1)

    return times, fluxs, full_ood


def get_4029_manual_mask(times, fluxs, cadencenos, get_full=0):

    # made using glue, with the manual_masking_20230617.glu session
    maskpaths = np.sort(
        glob(join(RESULTSDIR, 'tic4029_segments', '*manual_mask*.csv'))
    )

    # out of dip times, fluxs, cadenceno's
    ood_df = pd.concat((pd.read_csv(m) for m in maskpaths))
    ood_df = ood_df.sort_values(by='cadenceno')

    # full out of dip mask
    full_ood = np.in1d(cadencenos, np.array(ood_df.cadenceno))

    # sectors 18 and 19 out of dip
    s18s19_ood = full_ood & (times > 1750) & (times < 1900)

    # sectors 25 and 26 out of dip
    s25s26_ood = full_ood & (times > 1950) & (times < 2100)

    # sectors 53, 58, 59 out of dip
    s53s58s59_ood =  full_ood & (times > 2600)

    if get_full:
        return full_ood

    if get_split:
        return s18s19_ood, s25s26_ood, s53s58s59_ood, full_ood



def _get_lcpaths_given_ticid(ticid, lcpipeline, require_lc=1, cachedir=None):
    # ticid like "289840928"

    assert isinstance(ticid, str)
    if ticid.startswith("TIC"):
        ticid = ticid.lstrip("TIC_")
        ticid = ticid.lstrip("TIC ")
        ticid = ticid.lstrip("TIC")

    if ticid.startswith("TIC"):
        ticid_str = ticid
    else:
        ticid_str = f"TIC {ticid}"

    if lcpipeline not in 'tars':

        lcset = lk.search_lightcurve(ticid_str)

        if lcpipeline == 'spoc2min':
            sel = (lcset.author=='SPOC') & (lcset.exptime.value == 120)
        elif lcpipeline == 'qlp':
            sel = (lcset.author=='QLP')
        elif lcpipeline == 'cdips':
            sel = (lcset.author=='CDIPS')
        elif lcpipeline == 'tess-spoc':
            sel = (lcset.author=='TESS-SPOC')
        elif lcpipeline == 'spoc2min_tars':
            # TARS to be appended
            sel = (lcset.author=='SPOC') & (lcset.exptime.value == 120)
        else:
            raise NotImplementedError('allows: spoc2min, qlp, cdips, tess-spoc, spoc2min_tars')

        lcc = lcset[sel].download_all(download_dir=cachedir)

    if cachedir is None:
        cachedir = LKCACHEDIR

    if 'tars' in lcpipeline:
        from tehsors.tess import (
            get_tars_lcs_add_wrapunpopular, get_tars_lcs
        )
        tarslcdir = join(TARSCACHEDIR, f'tic_{ticid}')
        ADD_UNPOPULAR = 0
        if ADD_UNPOPULAR:
            time_dict, flux_dict = get_tars_lcs_add_wrapunpopular(ticid, tars_cache_dir=tarslcdir)
            ffi_dir = os.path.expanduser("~/.unpopular_cache")
            wrapunpoppaths = np.sort(glob(join(ffi_dir, f'TIC{ticid}*llc.csv')))
        else:
            time_dict, flux_dict = get_tars_lcs(ticid, cache_dir=tarslcdir)
        tarspaths = np.sort(glob(join(tarslcdir, f'tic_{ticid}_s*.csv')))

    if lcpipeline == 'spoc2min':
        lcpaths = glob(
            join(cachedir, f'tess*{ticid}*-s', f'tess*{ticid}*-s_lc.fits')
        )
    elif lcpipeline in ['qlp', 'cdips', 'tess-spoc']:
        lcpaths = glob(
            join(cachedir.replace("TESS","HLSP"), f'hlsp_{lcpipeline}*{ticid}*', f'*{ticid}*.fits')
        )
    elif lcpipeline == 'tars':
        if ADD_UNPOPULAR:
            lcpaths = list(tarspaths) + list(wrapunpoppaths)
        else:
            lcpaths = list(tarspaths)
    elif lcpipeline == 'spoc2min_tars':
        lcpaths = list(glob(
            join(cachedir, f'tess*{ticid}*-s', f'tess*{ticid}*-s_lc.fits')
        )) + list(tarspaths)

    if require_lc:
        msg = f'{ticid}: did not get LC'
        assert len(lcpaths) > 0, msg

    return lcpaths


def get_qlp_lcpaths(ticid):

    assert isinstance(ticid, str)

    from complexrotators.paths import QLPDIR
    # this CSV file has only the <500pc stars, which cuts time ~5x relative to
    # the full query
    lc_cache_path = join(
        QLPDIR, "s0001_to_s0055_QLP_ticid_path_merged_parallax_gt_2mas.csv"
    )

    delim = ","

    # use grep rather than reading the whole thing with pandas (much faster)
    colnames = os.popen( f'head -n1 {lc_cache_path}' ).read()
    colnames = colnames.rstrip('\n').split(delim)

    rowentry = os.popen( f'grep {ticid} {lc_cache_path}' ).read()

    if len(rowentry) >= 1:

        rowentry = rowentry.split("\n")[:-1]

        df = pd.DataFrame({})
        for ix, k in enumerate(colnames):
            df[k] = [r.split(",")[ix] for r in rowentry]

        lcpaths = list(df.path)
        lcpaths = [join(QLPDIR, l) for l in lcpaths]

        return lcpaths

    else:

        return []


def fix_cdips_lcpaths(lcpaths):
    """
    some look like
    /ar1/TESS/CDIPS/CDIPS/s0027/cam3_ccd3/hlsp_cdips_tess_ffi_gaiatwo0005196433038252282368-0027-cam3-ccd3_tess_v01_llc.fits
    perhaps a mast error?

    should be "s0027".  This fixes that.
    """

    basedirs = [os.path.dirname(l) for l in lcpaths]
    lcnames =  [os.path.basename(l) for l in lcpaths]

    outlcpaths = []

    # a weird name error exists in the cache csv file
    for dn,l in zip(basedirs, lcnames):
        a,b,c,d = l.split('-')
        if not b.startswith("s"):
            b = "s"+b
        outname = f"{a}-{b}-{c}-{d}"
        outpath = os.path.join(dn, outname)
        outlcpaths.append(outpath)

    return outlcpaths



def get_cdips_lcpaths(dr2_source_id):

    assert isinstance(dr2_source_id, str)

    from complexrotators.paths import CDIPSDIR
    # this CSV file has only the <500pc stars, which cuts time ~5x relative to
    # the full query
    lc_cache_path = join(
        CDIPSDIR, "hlsp_cdips_tess_ffi_s0001-s0055_tess_v01_catalog.csv"
    )

    delim = ";"

    # use grep rather than reading the whole thing with pandas (much faster)
    colnames = os.popen( f'head -n1 {lc_cache_path}' ).read()
    colnames = colnames.rstrip('\n').split(delim)

    rowentry = os.popen( f'grep {dr2_source_id} {lc_cache_path}' ).read()

    if len(rowentry) >= 1:

        rowentry = rowentry.split("\n")[:-1]

        df = pd.DataFrame({})
        for ix, k in enumerate(colnames):
            df[k] = [r.split(delim)[ix] for r in rowentry]

        lcpaths = list(df.name)
        lcpaths = [join(CDIPSDIR, "CDIPS", l) for l in lcpaths]

        lcpaths = fix_cdips_lcpaths(lcpaths)

        return lcpaths

    else:

        return []


def _get_local_lcpaths_given_ticid(ticid, lcpipeline, dr2_source_id=None):

    if lcpipeline == 'spoc2min':
        from complexrotators.paths import SPOCDIR
        lcpaths = glob(join(SPOCDIR, "sector-*", f"*{ticid}*.fits"))

    elif lcpipeline == 'qlp':
        lcpaths = get_qlp_lcpaths(ticid)

    elif lcpipeline == 'cdips':
        lcpaths = get_cdips_lcpaths(dr2_source_id)

    if len(lcpaths) == 0:
        print(f"Failed to find {lcpipeline} light curves for {ticid}")

    return lcpaths


def get_2min_cadence_spoc_tess_lightcurve(
    ticstr: str,
    ) -> list:
    """
    Args:
        ticstr: e.g., 'TIC 441420236'.  will be passed to lightkurve and used
        in naming plots.

    Returns:
        list containing individual sector-by-sector light curves.  If none are
        found, an empty list is returned.
    """

    # get the light curve
    lcset = lk.search_lightcurve(ticstr)

    if len(lcset) == 0:
        return []

    sel = (lcset.author=='SPOC') & (lcset.exptime.value == 120)
    lcc = lcset[sel].download_all()

    if lcc is None:
        return []

    # select only the two-minute cadence SPOC-reduced data; convert to a list.
    # note that this conversion approach works for any LightCurveCollection
    # returned by lightkurve -- no need to hand-pick the right ones.  the exact
    # condition below says "if the interval is between 119 and 121 seconds,
    # take it".
    lc_list = [_l for _l in lcc
         if
         _l.meta['ORIGIN']=='NASA/Ames'
         and
         np.isclose(
             120,
             np.nanmedian(np.diff(_l.remove_outliers().time.value))*24*60*60,
             atol=5
         )
    ]

    tstarts = [_lc.meta['TSTART'] for _lc in lcc]

    sorted_inds = np.argsort(tstarts)

    lc_list = [lc_list[i] for i in sorted_inds]

    return lc_list


def get_20sec_cadence_spoc_tess_lightcurve(
    ticstr: str,
    ) -> list:
    """
    Args:
        ticstr: e.g., 'TIC 441420236'.  will be passed to lightkurve and used
        in naming plots.

    Returns:
        list containing individual sector-by-sector light curves.  If none are
        found, an empty list is returned.
    """

    # get the light curve
    lcset = lk.search_lightcurve(ticstr)

    if len(lcset) == 0:
        return []

    sel = (lcset.author=='SPOC') & (lcset.exptime.value == 20)
    lcc = lcset[sel].download_all()

    if lcc is None:
        return []

    # select only the two-minute cadence SPOC-reduced data; convert to a list.
    # note that this conversion approach works for any LightCurveCollection
    # returned by lightkurve -- no need to hand-pick the right ones.  the exact
    # condition below says "if the interval is between 119 and 121 seconds,
    # take it".
    lc_list = [_l for _l in lcc
         if
         _l.meta['ORIGIN']=='NASA/Ames'
         and
         np.isclose(
             20,
             np.nanmedian(np.diff(_l.remove_outliers().time.value))*24*60*60,
             atol=5
         )
    ]

    return lc_list


def get_tic8_row(ticid, cachedir):

    import time

    cachepath = join(cachedir, f"TIC_{ticid}_mast_tic8_query.csv")

    if os.path.exists(cachepath):
        print(f'Found {cachepath}, returning.')
        return pd.read_csv(cachepath)

    ticstr = f"TIC {ticid}"
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
            print(f'TIC {ticid} failed initial MAST query with {e}, retrying {ix}/{MAX_ITER}...')
    if tic_data is None:
        raise ResolverError(f'TIC {ticid} failed to get MAST query.')

    t8_row = pd.DataFrame(tic_data.to_pandas().iloc[0]).T

    t8_row = t8_row.rename({c:f"tic8_{c}" for c in t8_row.columns},
                           axis='columns')

    t8_row.to_csv(cachepath, index=False)
    print(f'Cached {cachepath}')

    return t8_row


def fit_pseudocontinuum(
    wav, flx, lambda0,
    dlambda_cont=50.,
    exclude_wav_bands=None,
    poly_deg=2,
    n_iter=4,
    clip_hi=2.0,
    clip_lo=3.0
):
    """
    Fit a pseudo-continuum around a spectral line from many line-free pixels,
    rather than sampling a single anchor pixel.

    The fit is an iteratively sigma-clipped polynomial over the pixels with
    |wav - lambda0| < dlambda_cont, excluding `exclude_wav_bands`.  Clipping is
    asymmetric: positive outliers (emission components, cosmic rays) are
    rejected harder (clip_hi sigma) than negative ones (clip_lo sigma), so a
    time-variable emission feature cannot inflate the continuum estimate.

    Args:
        wav, flx: full-order wavelength [angstrom] and flux arrays (i.e. NOT
            pre-cut to the narrow display window around the line).
        lambda0: line rest wavelength [angstrom].
        dlambda_cont: half-width [angstrom] of the fitting window.
        exclude_wav_bands: list of (wav_lo, wav_hi) intervals [angstrom] to
            exclude from the fit (line core+wings, known blends).
        poly_deg: polynomial degree (2 allows for residual blaze curvature).
        n_iter: number of sigma-clipping iterations.
        clip_hi / clip_lo: positive / negative clip thresholds, in units of the
            robust (MAD-based) residual scatter.

    Returns dict with keys:
        coeffs: np.polyfit coefficients in x = (wav - lambda0)/dlambda_cont.
        poly: callable, evaluates the continuum at wavelength(s) [angstrom].
        cont_lambda0: continuum evaluated at lambda0.
        cont_frac_err: fractional statistical uncertainty of cont_lambda0
            (scatter of unclipped pixels about the fit / sqrt(N_kept)).
        n_kept: number of pixels in the final fit.
    """
    if exclude_wav_bands is None:
        exclude_wav_bands = []

    sel = np.abs(wav - lambda0) < dlambda_cont
    w, f = wav[sel], flx[sel]

    keep0 = np.isfinite(f) & np.isfinite(w)
    for (w0, w1) in exclude_wav_bands:
        keep0 &= ~((w >= w0) & (w <= w1))

    assert keep0.sum() > 5 * (poly_deg + 1), (
        'too few line-free pixels for pseudo-continuum fit'
    )

    # center/scale the abscissa for numerical conditioning
    x = (w - lambda0) / dlambda_cont

    keep = keep0.copy()
    for _ in range(n_iter):
        coeffs = np.polyfit(x[keep], f[keep], poly_deg)
        resid = f - np.polyval(coeffs, x)
        _r = resid[keep]
        sigma = 1.4826 * np.nanmedian(np.abs(_r - np.nanmedian(_r)))
        newkeep = keep0 & (resid < clip_hi * sigma) & (resid > -clip_lo * sigma)
        if newkeep.sum() < max(2 * (poly_deg + 1), 0.3 * keep0.sum()):
            break
        if (newkeep == keep).all():
            keep = newkeep
            break
        keep = newkeep

    coeffs = np.polyfit(x[keep], f[keep], poly_deg)
    resid = f[keep] - np.polyval(coeffs, x[keep])
    cont_lambda0 = float(np.polyval(coeffs, 0.))
    cont_frac_err = float(
        np.std(resid) / np.sqrt(keep.sum()) / np.abs(cont_lambda0)
    )

    poly = lambda _w: np.polyval(coeffs, (_w - lambda0) / dlambda_cont)

    return {
        'coeffs': coeffs,
        'poly': poly,
        'cont_lambda0': cont_lambda0,
        'cont_frac_err': cont_frac_err,
        'n_kept': int(keep.sum()),
    }


# Wavelength intervals excluded from the pseudo-continuum fit, per line.  The
# line core+wings exclusion (CONT_CORE_MASK_KMS, symmetric in velocity so that
# red- and blue-shifted emission are treated identically) is applied on top of
# these.  Vetted against the MagE night-median spectra for Hα (Ca I 6572.78
# falls inside the core+wings mask, but is listed explicitly in case the core
# mask is ever narrowed).  NOT yet vetted for Hβ.
CONT_CORE_MASK_KMS = 600.
CONT_EXTRA_MASKS = {
    'Hα': [(6572.78 - 1.5, 6572.78 + 1.5)],  # Ca I 6572.78
    'Hβ': [],
}


def get_specriver_data(
    ticid, linestr, specdir=os.path.expanduser('~/Dropbox/proj/cpv/data/spectra/HIRES/'),
    dlambda=20,
    usespectype='deblazed',
    gaussian_filter_sigma=2,
    datestr=None,
    normmethod=None,
    dlambda_cont=50.
):
    """
    Get flux vs wavelength cutouts for your desired line, given HIRES/MIKE/MagE
    .fits spectra in a directory.

    usespectype: 'deblazed' or 'reduced'

    normmethod: 'anchor', 'polycont', or None.
        'anchor': legacy behavior; the continuum is a single pixel of a
            sigma=5-pixel Gaussian-smoothed cutout at a fixed velocity offset
            (normatvel).  For MagE (23.4 km/s/pix) that kernel is 117 km/s
            wide, so nearby emission components leak into the normalization
            (see mage_continuum_normalization_assessment.txt).
        'polycont': fit a sigma-clipped polynomial pseudo-continuum over
            +/- dlambda_cont angstrom (line core+wings and known blends
            masked, symmetric in velocity), and rectify the cutout by the
            fitted polynomial.  norm_flx is the continuum at line center.
        None: resolves to 'polycont' for MagE and 'anchor' otherwise (the
            HIRES/MIKE kernels are 6.6/11 km/s, where the anchor pathology is
            small, and this keeps those published products reproducible).

    Returns:
        specpaths, spectimes, xvals, yvals, yvalsnonorm, norm_flxs,
        norm_flx_errs
        (norm_flx_errs: fractional continuum uncertainties for 'polycont',
         else NaN.)
    """

    from cdips_followup.spectools import read_hires
    from cdips_followup.spectools import read_mike
    from cdips_followup.spectools import read_mage

    # assumption: you have downloaded the "reduced" and the "deblazed" data
    # products from JUMP.  "blaze", generated by separate analysis scripts, is
    # optional, but needed if you want to fine-tune the blaze subtraction.

    if ticid in ['402980664', '141146667']:
        reducdir = join(specdir, f"TIC{ticid}_REDUCED")
        if ticid == '402980664':
            reducdir = join(specdir, f"TIC{ticid}_RDX")
        deblazedir = join(specdir, f"TIC{ticid}_DEBLAZED")
        blazedir = join(specdir, f"TIC{ticid}_BLAZE") # blaze function generated by get_tic1411_halpha_blaze.py here

    if usespectype == 'deblazed':
        datadir = deblazedir
    elif usespectype == 'reduced':
        datadir = reducdir
    elif usespectype == 'mike':
        datadir = os.path.expanduser('~/Dropbox/proj/cpv/data/spectra/MIKE/RAW/20250102/Final-Science-Products')
    elif 'mage' in usespectype:
        datadir = os.path.expanduser('~/Dropbox/proj/cpv/data/spectra/MAGE_REDUX/Final-Products')
    assert os.path.exists(datadir)

    if ticid in ['402980664', '141146667']:
        timecsvpath = join(datadir, f'tic{ticid}_trv.csv') # NOTE: this must be downloaded from JUMP
        timedf = pd.read_csv(timecsvpath, comment='#')
    elif ticid == '300651846' and usespectype == 'mike':
        timecsvpath = join(datadir, 'tic3006_bjd_times.csv') # NOTE: made by compute_bjd_times.py
        timedf = pd.read_csv(timecsvpath)
    elif ticid == '300651846' and usespectype == 'mage':
        timecsvpath = join(datadir, 'tic3006_mage_bjd_times.csv') # NOTE: made by compute_bjd_times.py
        timedf = pd.read_csv(timecsvpath)

    normpointdict = {
        '141146667': 700,
        '402980664': 350,
        '300651846': 89.65*5.5,
    }
    datedict = {
        '141146667': 'j537',
        '402980664': 'j546',
        '300651846': None # For TIC 300651846, datestr must be passed
    }
    radecdict = {
        # ra,dec
        '141146667': [166.31, 59.25],
        '402980664': [16.733, 80.459],
        '300651846': [114.360833, -66.752500]
    }
    assert str(ticid) in datedict

    if datestr is None:
        datestr = datedict[str(ticid)]

    if ticid in ['402980664', '141146667']:
        sel = timedf.observation_id.str.contains(datestr)
        seltimedf = timedf[sel]
    elif ticid == '300651846' and usespectype == 'mike':
        seltimedf = timedf
    elif ticid == '300651846' and usespectype == 'mage':
        sel = timedf.observation_id.astype(str).str.contains(datestr)
        seltimedf = timedf[sel]
    else:
        raise NotImplementedError

    assert linestr in ['Hα', 'Hβ', 'Hγ', 'Hδ', 'He', 'Ca[K]', 'K', 'Li']
    if ticid in ['402980664', '141146667']: # HIRES
        infodict = {
            # λ0, chip, order, normpoint(km/s)
            'Hα': [6562.8, 'i', 0, normpointdict[ticid]],
            'Hγ': [4340.47, 'b', 15, None],
            'Hδ': [4101.75, 'b', 10, None],
            'Ca[K]': [3933.66, 'b', 6, None],
            'He': [5875.62, 'r', 10, None],
            'K': [7699, 'i', 8, None],
            'Li': [6708, 'i', 1, None],
        }
    elif ticid in ['300651846'] and usespectype == 'mike':
        infodict = {
            # λ0, chip, order, normpoint(km/s)
            'Hα': [6562.8, 'red', 15, normpointdict[ticid]],
            'Hβ': [4861.33, 'blue', 2, normpointdict[ticid]],
            'Hγ': [4340.47, 'b', 15, None],
            'Hδ': [4101.75, 'b', 10, None],
        }
    elif ticid in ['300651846'] and usespectype == 'mage':
        infodict = {
            # λ0, chip, order, normpoint(km/s)
            'Hα': [6562.8, None, 3, normpointdict[ticid]],
            'Hβ': [4861.33, None, 7, normpointdict[ticid]],
            'Hγ': [4340.47, None, 8, None],
            'Hδ': [4101.75, None, 9, None],
        }

    if linestr != 'Hα' and usespectype == 'reduced':
        print(
            42*'!'+'\n'+
            'WARNING: BLAZE FUNCTION NOT REMOVED\n'+
            42*'!'+'\n'+
            42*'!'+'\n'+
            'WARNING: BLAZE FUNCTION NOT REMOVED\n'+
            42*'!'+'\n'
        )

    λ0 = infodict[linestr][0]
    λmin, λmax = λ0 - dlambda, λ0 + dlambda
    chip = infodict[linestr][1]
    order = infodict[linestr][2]
    normatvel = infodict[linestr][3]

    if normmethod is None:
        normmethod = 'polycont' if usespectype == 'mage' else 'anchor'
    assert normmethod in ['anchor', 'polycont']
    print(f'get_specriver_data: normmethod={normmethod}')

    if normmethod == 'polycont':
        _c_kms = const.c.to(u.km/u.s).value
        _dw_core = λ0 * CONT_CORE_MASK_KMS / _c_kms
        cont_exclude_bands = (
            [(λ0 - _dw_core, λ0 + _dw_core)]
            + CONT_EXTRA_MASKS.get(linestr, [])
        )

    def get_vel(wav, wav0):
        deltawvlen = ( wav - wav0 )
        delta_v = const.c * (deltawvlen / wav0)
        delta_v_kms = delta_v.to(u.km/u.s)
        return delta_v_kms.value

    if ticid in ['402980664', '141146667']: # HIRES
        specpaths = np.sort(glob(join(datadir, f"{chip}{datestr}*.fits")))
    elif ticid in ['300651846'] and usespectype == 'mike':
        specpaths = np.sort(glob(join(datadir, f"tic3006{chip}*.fits")))
    elif ticid in ['300651846'] and usespectype == 'mage':
        # NB: the historical exclusion of the first 20250119 exposure
        # (spec1184) is now applied INSIDE the loop, after the stale-duplicate
        # check, so that spec1185 (a stale copy of spec1184's data) is still
        # recognized and dropped.
        if datestr == '20250118':
            specpaths = np.sort(glob(join(datadir, f"tic3006_spec0*.fits")))
        elif datestr == '20250119':
            specpaths = np.sort(glob(join(datadir, f"tic3006_spec1*.fits")))
        elif datestr == '20250120':
            specpaths = np.sort(glob(join(datadir, f"tic3006_spec2*.fits")))
        elif datestr == '20250121':
            specpaths = np.sort(glob(join(datadir, f"tic3006_spec3*.fits")))
    assert len(specpaths) > 0

    yvals, xvals, spectimes = [], [], []
    yvalsnonorm = []
    norm_flxs = []
    norm_flx_errs = []
    keptpaths = []

    # For MagE, match each spectrum's time by filename (the CSV's row order
    # need not track the kept-file order once duplicates are dropped).
    if usespectype == 'mage':
        _mage_time_by_base = {
            os.path.basename(f): t for f, t in
            zip(timedf['file'], timedf['bjd_tdb_midtime'])
        }
        _prev_order_flx = None
        _n_dropped_dupes = 0

    if ticid in ['402980664', '141146667']: # HIRES
        # NOTE hack in case you opt for "REDUCED" only..
        deblazepaths = np.sort(glob(join(deblazedir, f"{chip}{datestr}*.fits")))
        reducpaths = np.sort(glob(join(reducdir, f"{chip}{datestr}*.fits")))

        if len(deblazepaths) == 0 and ticid == '402980664':
            deblazepaths = reducpaths
        defaultdeblazepath = deblazepaths[0]
        _, wav_2d = read_hires(defaultdeblazepath, is_registered=0, return_err=0)
    else:
        reducpaths = specpaths

    for _ix, (specpath, reducpath) in enumerate(zip(specpaths, reducpaths)):
        print(_ix, specpath)

        obs_id = os.path.basename(specpath).rstrip(".fits").lstrip(chip)

        if ticid in ['402980664', '141146667']: # HIRES
            try:
                flx_2d, wav_2d = read_hires(specpath, is_registered=0, return_err=0)
            except IndexError:
                # odd edge case for j537.174 blue chip...
                print('caught index error and pulling just flux...')
                hdul = fits.open(specpath)
                flx_2d = hdul[0].data
                hdul.close()
        elif ticid in ['300651846'] and usespectype == 'mike':
            flx_2d, wav_2d = read_mike(specpath)
        elif ticid in ['300651846'] and usespectype == 'mage':
            flx_2d, wav_2d = read_mage(specpath)

        start = 10
        end = -10
        flx, wav = flx_2d[order, start:end], wav_2d[order, start:end]

        # MAGE_REDUX contains exposures whose extracted spectra are stale
        # numerical copies of the preceding exposure (headers/times differ;
        # data identical to ~1e-7 relative, vs >~8% for genuinely distinct
        # neighbors).  Drop them: they would be plotted at phases at which
        # they were not taken, and double-weight the percentile template.
        if usespectype == 'mage':
            if (
                _prev_order_flx is not None
                and len(_prev_order_flx) == len(flx)
            ):
                _reldiff = (
                    np.nanmax(np.abs(flx - _prev_order_flx))
                    / np.nanmedian(np.abs(flx))
                )
                if _reldiff < 1e-4:
                    print(
                        f'!!! dropping {os.path.basename(specpath)}: '
                        f'order-{order} flux identical to previous kept '
                        f'exposure (reldiff {_reldiff:.1e}); stale duplicate '
                        'in MAGE_REDUX.'
                    )
                    _n_dropped_dupes += 1
                    continue
            _prev_order_flx = flx.copy()

            # historical exclusion of the first 20250119 exposure (was a
            # "[1:]" on the glob); applied after the duplicate check so that
            # spec1185 is still compared against spec1184.
            if (
                datestr == '20250119'
                and os.path.basename(specpath).startswith('tic3006_spec1184')
            ):
                print('!!! dropping tic3006_spec1184 (historical exclusion)')
                continue

        if (
            linestr == 'Hα' and usespectype == 'reduced' and
            ticid not in ['141146667', '402980664']
        ):
            if usespectype != 'mike':
                raise NotImplementedError(
                    'need to genearlize the blaze function approach, '
                    'or make it a kwarg'
                )

        if linestr == 'Hα' and ticid == '141146667':
            csvpath = join(
                blazedir, 'j537_ichip_order00_preferred_median_blaze_flux.csv'
            )
            badblaze = 'ij537.169'
            if usespectype == 'reduced':
                df = pd.read_csv(csvpath)
                blaze_flx = np.array(df.blz_flx)
                assert len(blaze_flx) == len(flx)
                #blz_flx_2d = reduc_flx_2d / deblz_flx_2d
                # -> deblz_flx_2d = reduc_flx_2d / blz_flx_2d
                flx = flx / blaze_flx
                print('fix deblaze')
            if usespectype == 'deblazed' and badblaze in specpath:
                hdul = fits.open(reducpath)
                flx_2d = hdul[0].data
                hdul.close()
                start = 10
                end = -10
                flx, wav = flx_2d[order, start:end], wav_2d[order, start:end]
                df = pd.read_csv(csvpath)
                blaze_flx = np.array(df.blz_flx)
                assert len(blaze_flx) == len(flx)
                # NOTE:
                # blz_flx_2d = reduc_flx_2d / deblz_flx_2d
                #  -> deblz_flx_2d = reduc_flx_2d / blz_flx_2d
                flx = flx / blaze_flx

        # full-order arrays, retained for the pseudo-continuum fit (the fit
        # window is wider than the display cutout).
        wav_order, flx_order = wav, flx

        sel = (wav > λmin) & (wav < λmax)
        wav = wav[sel]
        flx = flx[sel]

        vels = get_vel(wav, λ0)

        if normmethod == 'polycont':
            # Fit a sigma-clipped pseudo-continuum over +/- dlambda_cont with
            # the line core+wings and known blends masked, then rectify the
            # cutout by the fitted polynomial (this also removes the local
            # blaze/absorption slope).  norm_flx is the continuum at line
            # center, kept as the photon-count proxy used downstream.
            contfit = fit_pseudocontinuum(
                wav_order, flx_order, λ0,
                dlambda_cont=dlambda_cont,
                exclude_wav_bands=cont_exclude_bands
            )
            norm_flx = contfit['cont_lambda0']
            norm_flx_err = contfit['cont_frac_err']
            normed_flx = flx / contfit['poly'](wav)
        else:
            # Legacy 'anchor' normalization: flux at a single fixed velocity
            # offset.  This might mean the median flux in the selected
            # wavelength region.  It might mean a particular specified
            # velocity, if `normatvel` is non-null.  NB: there is no actual
            # smoothing happening in this step.
            if normatvel is None:
                norm_flx = np.nanmedian(flx)
            else:
                fn = lambda x: gaussian_filter1d(x, sigma=5)
                norm_flx = fn(flx)[ np.argmin(abs(vels - normatvel)) ]
            norm_flx_err = np.nan
            normed_flx = flx / norm_flx

        # NOTE : by default you are doing some gaussian smoothing...
        if gaussian_filter_sigma is not None:
            fn = lambda x: gaussian_filter1d(
                x, sigma=gaussian_filter_sigma
            )
            yvals.append(fn(normed_flx))
            yvalsnonorm.append(fn(flx))
        else:
            yvals.append(normed_flx)
            yvalsnonorm.append(flx)

        norm_flxs.append(norm_flx) # floating normalization due to seeing
        norm_flx_errs.append(norm_flx_err)
        xvals.append(vels)
        keptpaths.append(specpath)

        if ticid in ['402980664', '141146667']: # HIRES
            hl = fits.open(specpath)
            mjd = hl[0].header['MJD']
            hl.close()
            t = Time(mjd, format='mjd', scale='utc')

            ra = radecdict[str(ticid)][0]*u.deg
            dec = radecdict[str(ticid)][1]*u.deg
            t_bjd_tdb, bary_corr = astropy_utc_time_to_bjd_tdb(
                t, ra, dec, observatory='earthcenter', get_barycorr=1
            )
            print(f"barycorr is {bary_corr*24*60:.1f} minutes")

            # The above t_bjd_tdb is from the "MJD" header keyword... I am not even
            # sure if this is file creation time, midtime, or what.  Use the jump
            # time instead, which I think is photon weighted

            # spectimes.append(t_bjd_tdb - 2457000) # to TJD

            jumptime_bjd_tdb = seltimedf.loc[
                seltimedf.observation_id==obs_id, 'bjd'
            ].iloc[0]
            spectimes.append(jumptime_bjd_tdb - 2457000) # to TJD

        elif ticid == '300651846' and usespectype == 'mage':
            # match by filename, so dropped duplicates cannot misalign times
            _base = os.path.basename(specpath)
            assert _base in _mage_time_by_base, (
                f'{_base} missing from tic3006_mage_bjd_times.csv'
            )
            spectimes.append(_mage_time_by_base[_base] - 2457000)

        elif ticid == '300651846':
            spectimes = seltimedf.bjd_tdb_midtime - 2457000

    if usespectype == 'mage' and _n_dropped_dupes > 0:
        print(
            f'!!! get_specriver_data: dropped {_n_dropped_dupes} stale '
            f'duplicate exposure(s) for {datestr}; '
            f'{len(keptpaths)} independent spectra remain.'
        )

    specpaths = np.array(keptpaths)

    # ragged array handling bc of variable wavelength solutions
    from copy import deepcopy
    xvals_ucount = np.unique([len(i) for i in xvals])
    if len(xvals_ucount) > 1:
        assert np.max(np.diff([len(i) for i in xvals])) == 1
        big = max(xvals_ucount)
        smol = min(xvals_ucount)
        _xvals, _yvals = [], []
        for _x, _y in zip(xvals, yvals):
            if len(_x) == smol:
                _xvals.append(_x)
                _yvals.append(_y)
            else:
                _xvals.append(_x[:-1])
                _yvals.append(_y[:-1])
        xvals = deepcopy(_xvals)
        yvals = deepcopy(_yvals)

    return (specpaths, np.array(spectimes), xvals, yvals, yvalsnonorm,
            norm_flxs, norm_flx_errs)


def get_bochanski2007_m_standard(
    sptype=None, activestr='nactive'
    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Bochanski+2007 constructed M dwarf spectral templates from SDSS DR3.
    They span M0V to M9V (really even L0).

    They include both "active" and "non-active" stars, classified based on
    amount of H-alpha.  Optionally, a composite of both (highest SN)

    kwargs:

        sptype (str): e.g. "M2".
        activestr (str): "active", "nactive", or "all"

    returns:

        tuple of wavelength and flux from the spectral standard.
    """

    def build_wavelength_array(header: fits.Header, n_pixels: int) -> np.ndarray:
        """Build a wavelength array from linear WCS keywords in a FITS header.

        Args:
            header (fits.Header): FITS header containing CRVAL1, CRPIX1, and CD1_1.
            n_pixels (int): Number of pixels along the wavelength axis.

        Returns:
            np.ndarray: 1D array of wavelengths corresponding to each pixel.
        """
        crval1 = header["CRVAL1"]
        crpix1 = header["CRPIX1"]
        cd1_1 = header["CD1_1"]

        pixel_indices = np.arange(n_pixels, dtype=float)
        # Note: FITS pixels are 1-based, so we do (pixel + 1 - crpix1)
        wavelength = crval1 + (pixel_indices + 1.0 - crpix1) * cd1_1

        return wavelength


    def read_spectrum(filename: str) -> Tuple[np.ndarray, np.ndarray]:
        """Read a multi-plane 1D spectrum from a FITS file and build the wavelength.

        This assumes that the first axis is wavelength and the second axis
        contains the different flux planes (e.g., mean flux, median flux, S/N, etc).

        Args:
            filename (str): Path to the FITS file.

        Returns:
            Tuple[np.ndarray, np.ndarray]:
                A tuple of (wavelength, flux_data) where:
                - wavelength is a 1D array of wavelengths.
                - flux_data is a 2D array with shape (n_planes, n_pixels).
        """
        with fits.open(filename) as hdul:
            header = hdul[0].header
            flux_data = hdul[0].data

        n_pixels = flux_data.shape[1]
        wavelength = build_wavelength_array(header, n_pixels)

        # return median-normalized coadded flux
        return wavelength, flux_data[1,:]

    assert sptype in [f"M{i}" for i in range(10)], f"Invalid spectral type: {sptype}"

    datadir = join(DATADIR, 'spectra', 'SDSS', 'Bochanski2007')

    namestr = f'{sptype.lower()}.{activestr}.na.k.fits'
    fitspath = join(datadir, namestr)
    assert os.path.exists(fitspath)

    wav, flx = read_spectrum(fitspath)

    return wav, flx
