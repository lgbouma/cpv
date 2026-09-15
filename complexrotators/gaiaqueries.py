"""
Query Gaia through a TAP mirror, as a drop-in replacement for
cdips.utils.gaiaqueries.given_source_ids_get_gaia_data.

Motivation: the ESA archive (gea.esac.esa.int) intermittently fails on its
/tap-server/tap/sync endpoint -- queries either hang or return HTTP 500, while
/tap-server/tap/availability keeps reporting 200.  When that happens,
plot_cpvvetter silently bails out before writing its PDF.

The DR2 tables are static, so any mirror serves identical values.  The ARI
Heidelberg mirror carries both gaiadr2.gaia_source and gaiadr2.ruwe.

Contents:
    | given_source_ids_get_gaia_data_mirror
    | given_source_ids_get_gaia_data_fallback
"""
#############
## LOGGING ##
#############
import logging
from complexrotators import log_sub, log_fmt, log_date_fmt

LOGGER = logging.getLogger(__name__)
LOGDEBUG = LOGGER.debug
LOGINFO = LOGGER.info
LOGWARNING = LOGGER.warning
LOGERROR = LOGGER.error

#############
## IMPORTS ##
#############
import os
from os.path import join
import numpy as np, pandas as pd

# ARI Heidelberg; an alternative is "https://gaia.aip.de/tap", but as of
# 2026/09 that one does not serve gaiadr2.ruwe.
GAIA_MIRROR_URL = "https://gaia.ari.uni-heidelberg.de/tap"

GAIACACHEDIR = join(os.path.expanduser("~"), ".gaia_cache")


def given_source_ids_get_gaia_data_mirror(
    source_ids, runid, n_max=5, overwrite=False,
    enforce_all_sourceids_viable=True, gaia_datarelease='gaiadr2',
    getdr2ruwe=False, mirror_url=GAIA_MIRROR_URL, cachedir=GAIACACHEDIR
    ):
    """
    Get Gaia rows for the given source_ids, via a TAP mirror.

    The signature mirrors the subset of
    cdips.utils.gaiaqueries.given_source_ids_get_gaia_data that
    plot_cpvvetter uses, so that it can be swapped in directly.

    Args:

        source_ids: np.ndarray of np.int64 Gaia source_ids.

        runid: string used to name the cache file.

        n_max: maximum number of rows to return.

        overwrite: if True, requery even when a cache file exists.

        enforce_all_sourceids_viable: if True, raise when the number of rows
            returned does not match the number of source_ids requested.

        gaia_datarelease: e.g. "gaiadr2".

        getdr2ruwe: if True, query the {gaia_datarelease}.ruwe table rather
            than {gaia_datarelease}.gaia_source.

    Returns:

        pd.DataFrame of the matching rows.
    """

    import pyvo

    if not os.path.exists(cachedir):
        os.makedirs(cachedir, exist_ok=True)

    tablename = 'ruwe' if getdr2ruwe else 'gaia_source'
    cachepath = join(
        cachedir, f"mirror_{runid}_{gaia_datarelease}_{tablename}.csv"
    )

    if os.path.exists(cachepath) and not overwrite:
        LOGINFO(f"Found {cachepath}")
        return pd.read_csv(cachepath)

    source_ids = np.atleast_1d(np.array(source_ids)).astype(np.int64)
    idstr = ", ".join([str(s) for s in source_ids])

    query = (
        f"SELECT TOP {n_max} * "
        f"FROM {gaia_datarelease}.{tablename} "
        f"WHERE source_id IN ({idstr})"
    )

    LOGINFO(f"Querying {mirror_url} for {runid} ({tablename})...")
    service = pyvo.dal.TAPService(mirror_url)
    result = service.search(query)

    df = result.to_table().to_pandas()

    if enforce_all_sourceids_viable:
        if len(df) != len(source_ids):
            errmsg = (
                f"{runid}: got {len(df)} rows from {mirror_url} for "
                f"{len(source_ids)} source_ids."
            )
            raise ValueError(errmsg)

    df.to_csv(cachepath, index=False)
    LOGINFO(f"Made {cachepath}")

    return df


def given_source_ids_get_gaia_data_fallback(
    source_ids, runid, n_max=5, overwrite=False,
    enforce_all_sourceids_viable=True, savstr='', which_columns='*',
    table_name='gaia_source', gaia_datarelease='gaiadr2', getdr2ruwe=False,
    mirror_first=True
    ):
    """
    Get Gaia rows, trying the mirror and the ESA archive in turn.

    Same call signature as cdips.utils.gaiaqueries.given_source_ids_get_gaia_data.

    mirror_first: if True (default), the TAP mirror is tried first and the ESA
        archive is the fallback.  The DR2 tables are static, so the two return
        identical values; the mirror is the default only because ESA has been
        unreliable, and responds in ~0.5s rather than timing out for minutes.
        Set False to restore ESA-first behaviour.
    """

    from cdips.utils.gaiaqueries import given_source_ids_get_gaia_data

    def _mirror():
        return given_source_ids_get_gaia_data_mirror(
            source_ids, runid, n_max=n_max, overwrite=overwrite,
            enforce_all_sourceids_viable=enforce_all_sourceids_viable,
            gaia_datarelease=gaia_datarelease, getdr2ruwe=getdr2ruwe
        )

    def _esa():
        return given_source_ids_get_gaia_data(
            source_ids, runid, n_max=n_max, overwrite=overwrite,
            enforce_all_sourceids_viable=enforce_all_sourceids_viable,
            savstr=savstr, which_columns=which_columns,
            table_name=table_name, gaia_datarelease=gaia_datarelease,
            getdr2ruwe=getdr2ruwe
        )

    if mirror_first:
        first, second, firstname, secondname = _mirror, _esa, 'mirror', 'ESA'
    else:
        first, second, firstname, secondname = _esa, _mirror, 'ESA', 'mirror'

    try:
        return first()
    except Exception as e:
        LOGWARNING(
            f"{runid}: {firstname} Gaia query failed ({e}); trying {secondname}."
        )

    return second()
