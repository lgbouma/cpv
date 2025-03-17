"""
Contents:

General:

    | get_tic8_row
    | get_2min_cadence_spoc_tess_lightcurve
    | get_20sec_cadence_spoc_tess_lightcurve

    | _get_lcpaths_given_ticid
    | _get_local_lcpaths_given_ticid
    | _get_lcpaths_fromlightkurve_given_ticid

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
from complexrotators.paths import LKCACHEDIR, LOCALDIR, RESULTSDIR, DATADIR

from astropy.io import fits
from astropy import units as u, constants as const
from astropy.time import Time
from scipy.ndimage import gaussian_filter1d

from astroquery.exceptions import ResolverError
from astroquery.mast import Catalogs

from cdips.utils.lcutils import astropy_utc_time_to_bjd_tdb

def get_cqv_search_sample():
    localdir = '/Users/luke/local/SPOCLC'
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



def _get_lcpaths_given_ticid(ticid):
    # TODO: this getter will need to be updated when running at scale on wh1
    SPOCDIR = "/Users/luke/local/SPOCLC"
    lcpaths = glob(join(SPOCDIR, "lightcurves", f"*{ticid}*.fits"))

    if len(lcpaths) == 0:
        p = subprocess.call([
            "scp", f"luke@wh1:/ar1/TESS/SPOCLC/sector*/*{ticid}*.fits",
            join(SPOCDIR, "lightcurves")
        ])
        assert 0
    return lcpaths


def _get_lcpaths_fromlightkurve_given_ticid(ticid, lcpipeline, require_lc=1,
                                            cachedir=None):
    # ticid like "289840928"

    assert isinstance(ticid, str)
    if ticid.startswith("TIC"):
        ticid_str = ticid
    else:
        ticid_str = f"TIC {ticid}"

    lcset = lk.search_lightcurve(ticid_str)

    if lcpipeline == 'spoc2min':
        sel = (lcset.author=='SPOC') & (lcset.exptime.value == 120)
    elif lcpipeline == 'qlp':
        sel = (lcset.author=='QLP')
    elif lcpipeline == 'cdips':
        sel = (lcset.author=='CDIPS')
    elif lcpipeline == 'tess-spoc':
        sel = (lcset.author=='TESS-SPOC')
    lcc = lcset[sel].download_all(download_dir=cachedir)

    if cachedir is None:
        cachedir = LKCACHEDIR

    if lcpipeline == 'spoc2min':
        lcpaths = glob(
            join(cachedir, f'tess*{ticid}*-s', f'tess*{ticid}*-s_lc.fits')
        )
    elif lcpipeline in ['qlp', 'cdips', 'tess-spoc']:
        lcpaths = glob(
            join(cachedir.replace("TESS","HLSP"), f'hlsp_{lcpipeline}*{ticid}*', f'*{ticid}*.fits')
        )

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


def get_specriver_data(
    ticid, linestr, specdir='/Users/luke/Dropbox/proj/cpv/data/spectra/HIRES/',
    dlambda=20,
    usespectype='deblazed',
    gaussian_filter_sigma=2
):
    """
    Get flux vs wavelength cutouts for your desired line, given HIRES .fits
    spectra in a directory.

    usespectype: 'deblazed' or 'reduced'
    """

    from cdips_followup.spectools import read_hires

    # assumption: you have downloaded the "reduced" and the "deblazed" data
    # products from JUMP.  "blaze", generated by separate analysis scripts, is
    # optional, but needed if you want to fine-tune the blaze subtraction.
    reducdir = join(specdir, f"TIC{ticid}_REDUCED")
    deblazedir = join(specdir, f"TIC{ticid}_DEBLAZED")
    blazedir = join(specdir, f"TIC{ticid}_BLAZE") # blaze function generated by get_tic1411_halpha_blaze.py here

    if usespectype == 'deblazed':
        datadir = deblazedir
    elif usespectype == 'reduced':
        datadir = reducdir
    assert os.path.exists(datadir)

    timecsvpath = join(datadir, f'tic{ticid}_trv.csv') # NOTE: this must be downloaded from JUMP
    timedf = pd.read_csv(timecsvpath, comment='#')

    normpointdict = {
        '141146667': 700,
        '402980664': 350
    }
    datedict = {
        '141146667': 'j537',
        '402980664': 'j531'
    }
    radecdict = {
        # ra,dec
        '141146667': [166.31, 59.25],
        '402980664': [16.733, 80.459]
    }
    assert str(ticid) in datedict
    datestr = datedict[str(ticid)]

    sel = timedf.observation_id.str.contains(datedict[ticid])
    seltimedf = timedf[sel]

    assert linestr in ['Hα', 'Hγ', 'Hδ', 'He', 'Ca[K]', 'K', 'Li']
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

    def get_vel(wav, wav0):
        deltawvlen = ( wav - wav0 )
        delta_v = const.c * (deltawvlen / wav0)
        delta_v_kms = delta_v.to(u.km/u.s)
        return delta_v_kms.value

    specpaths = np.sort(glob(join(datadir, f"{chip}{datestr}*.fits")))
    assert len(specpaths) > 0

    yvals, xvals, spectimes = [], [], []
    yvalsnonorm = []
    norm_flxs = []

    # NOTE hack in case you opt for "REDUCED" only..
    deblazepaths = np.sort(glob(join(deblazedir, f"{chip}{datestr}*.fits")))
    reducpaths = np.sort(glob(join(reducdir, f"{chip}{datestr}*.fits")))

    defaultdeblazepath = deblazepaths[0]
    _, wav_2d = read_hires(defaultdeblazepath, is_registered=0, return_err=0)

    for specpath, reducpath in zip(specpaths, reducpaths):

        obs_id = os.path.basename(specpath).rstrip(".fits").lstrip(chip)

        try:
            flx_2d, wav_2d = read_hires(specpath, is_registered=0, return_err=0)
        except IndexError:
            # odd edge case for j537.174 blue chip...
            print('caught index error and pulling just flux...')
            hdul = fits.open(specpath)
            flx_2d = hdul[0].data
            hdul.close()
        start = 10
        end = -10
        flx, wav = flx_2d[order, start:end], wav_2d[order, start:end]

        if linestr == 'Hα' and usespectype == 'reduced' and ticid != '141146667':
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

        sel = (wav > λmin) & (wav < λmax)
        wav = wav[sel]
        flx = flx[sel]

        vels = get_vel(wav, λ0)

        # Normalize flux to 1 at whatever normatvel is given at.
        # This might mean the median flux in the selected wavelength
        # region.  It might mean a particular specified velocity, if
        # `normatvel` is non-null.  NB: there is no actual smoothing
        # happening in this step.
        if normatvel is None:
            norm_flx = np.nanmedian(flx)
        else:
            fn = lambda x: gaussian_filter1d(x, sigma=5)
            norm_flx = fn(flx)[ np.argmin(abs(vels - normatvel)) ]

        # NOTE : by default you are doing some gaussian smoothing...
        if gaussian_filter_sigma is not None:
            fn = lambda x: gaussian_filter1d(
                x, sigma=gaussian_filter_sigma
            )
            yvals.append(fn(flx/norm_flx))
            yvalsnonorm.append(fn(flx))
        else:
            yvals.append(flx/norm_flx)
            yvalsnonorm.append(flx)

        norm_flxs.append(norm_flx) # floating normalization due to seeing
        xvals.append(vels)

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

    return specpaths, np.array(spectimes), xvals, yvals, yvalsnonorm, norm_flxs


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
