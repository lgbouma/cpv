"""
Contents:
| get_2min_cadence_spoc_tess_lightcurve
| get_20sec_cadence_spoc_tess_lightcurve
| _get_lcpaths_given_ticid
"""
import numpy as np
import lightkurve as lk
from os.path import join
from glob import glob

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
