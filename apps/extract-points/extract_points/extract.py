"""Orchestrate digitization: panels -> calibration -> markers -> data points.

Produces, for every flux panel, the per-band digitized (BJD-offset, relative
flux) points, plus a tidy long-form table across all panels.

Note on the IR bands: Tanimoto et al. plot the IR series (K_s/J/H) with an
arbitrary downward display offset for clarity. We keep that offset so the
reproduced panels match the originals; the optical (I_C) series is on the
plotted relative-flux scale (baseline ~1.0).
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
import pandas as pd

from .image_io import load_rgb, dark_mask
from .panels import detect_flux_boxes, FluxBox
from .calibrate import calibrate, Calibration
from .markers import (
    clean_mask, split_band_masks, detect_markers_matched, bin_cadence)
from .config import PANELS, CSV_EPOCHS

# extraction methods: how to build the marker mask (red_fill) and how to reduce
# a band's pixels to (x_px, y_px) points.
METHODS = {
    # Option 1: matched-filter per-marker detection on the black-marker mask
    "per_marker": dict(reduce=detect_markers_matched, red_fill=False),
    # Option 2: cadence-matched binning (median y)
    "cadence": dict(reduce=bin_cadence, red_fill=False),
    # Option 3: per-marker detection after restoring red-model-occluded interiors
    "red_aware": dict(reduce=detect_markers_matched, red_fill=True),
}


@dataclass
class BandData:
    band: str
    is_optical: bool       # True -> filled circle (I_C); False -> cross (IR)
    bjd: np.ndarray        # BJD - offset
    flux: np.ndarray       # relative flux (IR includes display offset)


@dataclass
class PanelResult:
    figure: str
    row: int
    col: int
    epoch: str
    offset: int
    in_csv: bool
    calib: Calibration
    box: FluxBox                 # detected flux-box spine pixels (absolute)
    xlim: tuple[float, float]    # (left, right) plot x-range in BJD - offset
    ylim: tuple[float, float]    # (bottom, top) plot y-range in relative flux
    bands: list[BandData] = field(default_factory=list)


def _pts_to_arrays(pts, calib: Calibration, box: FluxBox):
    """Convert box-local (x, y) marker pixels to (BJD-offset, flux).

    extract_bands works on the cropped panel and returns box-local pixel
    coordinates; the calibration tick positions are in absolute image pixels,
    so we add the box origin before applying the transform."""
    if not pts:
        return np.array([]), np.array([])
    arr = np.array(pts)
    order = np.argsort(arr[:, 0])
    arr = arr[order]
    bjd = calib.px_to_bjd(arr[:, 0] + box.x0)
    flux = calib.px_to_flux(arr[:, 1] + box.y0)
    return bjd, flux


def extract_panel(rgb: np.ndarray, dark: np.ndarray, box: FluxBox,
                  method: str = "red_aware") -> PanelResult:
    meta = PANELS[(box.figure, box.row, box.col)]
    calib = calibrate(dark, box)
    spec = METHODS[method]
    mm = clean_mask(rgb, box, red_fill=spec["red_fill"])
    reduce_band = spec["reduce"]
    band_masks = split_band_masks(mm, len(meta["bands"]))

    # The plot's true axis limits are the data values at the box spines:
    # x0/x1 are the left/right spines, y0/y1 the top spine / flux-box bottom.
    # px_to_bjd returns BJD - offset (the plotted x units).
    xlim = (float(calib.px_to_bjd(box.x0)), float(calib.px_to_bjd(box.x1)))
    ylim = (float(calib.px_to_flux(box.y1)), float(calib.px_to_flux(box.y0)))

    res = PanelResult(box.figure, box.row, box.col, meta["epoch"],
                      meta["offset"], meta["epoch"] in CSV_EPOCHS, calib,
                      box, xlim, ylim)

    # bands[0] = I_C (optical, upper, circle); bands[1] = IR (lower, cross)
    for bi, bmask in enumerate(band_masks):
        bjd, flux = _pts_to_arrays(reduce_band(bmask), calib, box)
        res.bands.append(BandData(meta["bands"][bi], bi == 0, bjd, flux))
    return res


def extract_figure(figure: str, method: str = "red_aware") -> list[PanelResult]:
    rgb = load_rgb(figure)
    dark = dark_mask(rgb)
    return [extract_panel(rgb, dark, b, method)
            for b in detect_flux_boxes(figure, rgb)]


def extract_all(method: str = "red_aware") -> list[PanelResult]:
    out = []
    for fig in ("f12", "f13"):
        out.extend(extract_figure(fig, method))
    return out


def to_dataframe(results: list[PanelResult]) -> pd.DataFrame:
    rows = []
    for r in results:
        for bd in r.bands:
            for bjd, flux in zip(bd.bjd, bd.flux):
                rows.append(dict(
                    figure=r.figure, panel_row=r.row, panel_col=r.col,
                    epoch=r.epoch, bjd_offset=r.offset, in_csv=r.in_csv,
                    band=bd.band, marker="circle" if bd.is_optical else "cross",
                    bjd_minus_offset=round(float(bjd), 6),
                    rel_flux=round(float(flux), 6),
                    bjd=round(float(bjd) + r.offset, 6),
                ))
    return pd.DataFrame(rows)


def panel_metadata_dataframe(results: list[PanelResult]) -> pd.DataFrame:
    """One row per panel: the original figure's axis limits (recovered from the
    box spines via the calibration), the box pixel bounds, calibration quality,
    and per-band point counts. Lets the replot match the originals' xlim/ylim
    and supports QA of the digitization."""
    rows = []
    for r in results:
        n_optical = r.bands[0].bjd.size if r.bands else 0
        n_ir = r.bands[1].bjd.size if len(r.bands) > 1 else 0
        rows.append(dict(
            figure=r.figure, panel_row=r.row, panel_col=r.col,
            epoch=r.epoch, bjd_offset=r.offset, in_csv=r.in_csv,
            bands="+".join(bd.band for bd in r.bands),
            # plotted x units are BJD - offset; add offset for absolute BJD
            xlim_lo=round(r.xlim[0], 6), xlim_hi=round(r.xlim[1], 6),
            ylim_lo=round(r.ylim[0], 6), ylim_hi=round(r.ylim[1], 6),
            x_r2=round(r.calib.x_r2, 6), y_r2=round(r.calib.y_r2, 6),
            n_xticks=len(r.calib.x_ticks_px), n_yticks=len(r.calib.y_ticks_px),
            box_x0=r.box.x0, box_x1=r.box.x1, box_y0=r.box.y0, box_y1=r.box.y1,
            n_optical=n_optical, n_ir=n_ir,
        ))
    return pd.DataFrame(rows)
