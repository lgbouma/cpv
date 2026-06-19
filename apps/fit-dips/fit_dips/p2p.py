"""
Vendored point-to-point RMS estimator.

Copied from cdips.utils.lcutils.p2p_rms (Bouma et al.) because the installed
`cdips` package is not importable in this environment.  See Section 2.6 of
Nardiello+2020 for the definition.
"""
import numpy as np


def p2p_rms(flux):
    """Point-to-point RMS: robust scatter insensitive to intrinsic variability.

    Computed from the 84th-16th percentile spread of the first differences
    delta F_j = F_j - F_{j+1}.
    """
    flux = np.asarray(flux, dtype=float)
    dflux = np.diff(flux)
    med_dflux = np.nanmedian(dflux)

    up_p2p = (
        np.nanpercentile(np.sort(dflux - med_dflux), 84)
        - np.nanpercentile(np.sort(dflux - med_dflux), 50)
    )
    lo_p2p = (
        np.nanpercentile(np.sort(dflux - med_dflux), 50)
        - np.nanpercentile(np.sort(dflux - med_dflux), 16)
    )

    return float(np.mean([up_p2p, lo_p2p]))
