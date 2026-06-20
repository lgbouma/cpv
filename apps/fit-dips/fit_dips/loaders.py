"""
Light-curve loaders for the heterogeneous CPV photometry formats.

Every loader returns (t, flux, flux_err) as finite numpy arrays with the flux
normalized to a median of 1.  flux_err is the instrument-provided error if any
(else None); note the fitting pipeline uses p2p_rms by default regardless.

Patterns mirror the proven readers in complexrotators/plotting.py and
drivers/plot_hires_lp12-502_mosaic.py.
"""
import os

import numpy as np
import pandas as pd

REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", ".."))


def _abspath(path):
    # Expand a leading "~" so data paths are portable across machines/home dirs
    # (paths are stored as e.g. "~/.lightkurve/..." rather than "/Users/<you>/").
    path = os.path.expanduser(path)
    if os.path.isabs(path):
        return path
    return os.path.join(REPO_ROOT, path)


def _finite_normalize(t, flux, flux_err=None, window=None):
    t = np.asarray(t, dtype=float)
    flux = np.asarray(flux, dtype=float)
    err = None if flux_err is None else np.asarray(flux_err, dtype=float)

    m = np.isfinite(t) & np.isfinite(flux)
    if err is not None:
        m &= np.isfinite(err)
    if window is not None:
        m &= (t >= window[0]) & (t <= window[1])
    t, flux = t[m], flux[m]
    if err is not None:
        err = err[m]

    order = np.argsort(t)
    t, flux = t[order], flux[order]
    if err is not None:
        err = err[order]

    med = np.nanmedian(flux)
    flux = flux / med
    if err is not None:
        err = err / med
    return t, flux, err


def load_muscat_csv(path, window=None, **_):
    df = pd.read_csv(_abspath(path))
    return _finite_normalize(df["BJD_TDB"], df["Flux"], df["Err"], window)


def load_tierras_csv(path, window=None, **_):
    df = pd.read_csv(_abspath(path))
    return _finite_normalize(
        df["BJD TDB"], df["Target Relative Flux"],
        df["Target Relative Flux Error"], window)


def load_keplercam_dat(path, window=None, **_):
    df = pd.read_csv(_abspath(path), sep=r"\s+")
    return _finite_normalize(
        df["BJD_TDB_B"], df["rel_flux_T1_n"], df["rel_flux_err_T1_n"], window)


def load_tess_fits(path, window=None, fluxkey="SAP_FLUX", **_):
    from astropy.io import fits
    with fits.open(_abspath(path)) as hdul:
        d = hdul[1].data
        bjd = np.array(d["TIME"]) + 2457000.0
        flux = np.array(d[fluxkey], dtype=float)
        qual = np.array(d["QUALITY"])
    bad = 2 | 4 | 8 | 128  # SafeMode, CoarsePoint, EarthPoint, ManualExclude
    good = (qual & bad) == 0
    bjd, flux = bjd[good], flux[good]
    return _finite_normalize(bjd, flux, None, window)


def load_muscat2_fits(path, band, window=None, **_):
    """MuSCAT2 achromatic FITS: one HDU per band (flux_g/r/i, flux_z_s)."""
    from astropy.io import fits
    fitskey = "flux_z_s" if band == "z" else f"flux_{band}"
    with fits.open(_abspath(path)) as hl:
        d = hl[fitskey].data
        t = np.array(d["time_bjd"], dtype=float)
        flux = np.array(d["flux"], dtype=float)
    return _finite_normalize(t, flux, None, window)


def load_fourstar_xls(path, window=None, flux_offset=0.76, tcorr_min=60.0,
                      trim_after=None, **_):
    """FourStar NB2.09 aperture photometry (tab-separated .xls).

    Replicates the reduction in drivers/plot_fourstar_vs_tess.py: a +tcorr_min
    time correction, a -flux_offset additive correction to rel_flux_C1, and an
    optional end-trim dropping points with t > trim_after.
    """
    df = pd.read_csv(_abspath(path), sep="\t")
    t = np.array(df["BJD_TDB"], dtype=float) + tcorr_min / (60.0 * 24.0)
    f = np.array(df["rel_flux_C1"], dtype=float) - flux_offset
    if trim_after is not None:
        keep = t <= trim_after
        t, f = t[keep], f[keep]
    return _finite_normalize(t, f, None, window)


def load_tess_pkl(path, window=None, time_offset=0.0, **_):
    """TESS light curve from a cpv_periodsearch pickle (d['times'], d['fluxs']).

    ``time_offset`` is added to the stored times (e.g. 2457000 to convert BTJD
    to BJD_TDB).
    """
    import pickle
    with open(_abspath(path), "rb") as fh:
        d = pickle.load(fh)
    t = np.array(d["times"], dtype=float) + time_offset
    return _finite_normalize(t, np.array(d["fluxs"]), None, window)


def load_generic_csv(path, time_col, flux_col, fluxerr_col=None, window=None,
                     sep=",", time_scale=1.0, flux_in_mag=False, **_):
    """Generic CSV. ``time_scale`` multiplies the time column (e.g. 1/24 to
    convert hours to days, or a rotation period to convert phase to days).

    ``flux_in_mag``: treat ``flux_col`` as a (differential) magnitude and
    convert to relative flux via 10**(-0.4*mag) before the median-1
    normalization, so a brightness dip (mag increase) becomes a flux decrement.
    """
    df = pd.read_csv(_abspath(path), sep=sep)
    t = np.asarray(df[time_col], dtype=float) * time_scale
    flux = np.asarray(df[flux_col], dtype=float)
    if flux_in_mag:
        flux = 10.0 ** (-0.4 * flux)
    err = df[fluxerr_col] if fluxerr_col else None
    return _finite_normalize(t, flux, err, window)


LOADERS = {
    "muscat_csv": load_muscat_csv,
    "muscat2_fits": load_muscat2_fits,
    "tierras_csv": load_tierras_csv,
    "keplercam_dat": load_keplercam_dat,
    "fourstar_xls": load_fourstar_xls,
    "tess_fits": load_tess_fits,
    "tess_pkl": load_tess_pkl,
    "generic_csv": load_generic_csv,
}


_NON_LOADER_KEYS = ("loader", "mask_indices", "fold_period", "fold_ref")


def load_lightcurve(data_spec):
    """Dispatch on data_spec['loader'] and return (t, flux, flux_err).

    Optional ``mask_indices``: integer indices (into the time-sorted, returned
    array) of points to drop before use, e.g. ``[0]`` to mask the first point.
    After masking, flux (and error) are re-normalized to a median of 1.

    Optional ``fold_period`` (days): phase-fold the times onto a single cycle
    via ``t -> t - round((t - fold_ref) / P) * P``, then re-sort. This maps
    every rotation cycle in the (typically widened) window onto the cycle
    containing ``fold_ref``, so that data from many cycles overlay in
    flux-vs-phase. ``fold_ref`` (days) defaults to the median time; choose it so
    the labeled dip/flare windows stay within +/- P/2 of it (i.e. do not wrap).
    Flux is unchanged by folding (so the median-1 normalization is preserved);
    only the time coordinate and sort order change.
    """
    loader = data_spec.get("loader")
    if loader not in LOADERS:
        raise ValueError(f"unknown loader: {loader!r}")
    kwargs = {k: v for k, v in data_spec.items()
              if k not in _NON_LOADER_KEYS}
    window = kwargs.pop("window", None)
    t, flux, err = LOADERS[loader](window=window, **kwargs)

    mask_indices = data_spec.get("mask_indices")
    if mask_indices:
        keep = np.ones(t.size, dtype=bool)
        keep[np.array(mask_indices, dtype=int)] = False
        t, flux = t[keep], flux[keep]
        err = err[keep] if err is not None else None
        med = np.nanmedian(flux)
        flux = flux / med
        if err is not None:
            err = err / med

    fold_period = data_spec.get("fold_period")
    if fold_period:
        P = float(fold_period)
        fold_ref = data_spec.get("fold_ref")
        fold_ref = float(np.median(t)) if fold_ref is None else float(fold_ref)
        t = t - np.round((t - fold_ref) / P) * P
        order = np.argsort(t)
        t, flux = t[order], flux[order]
        if err is not None:
            err = err[order]
    return t, flux, err
