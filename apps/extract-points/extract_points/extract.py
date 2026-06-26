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
from .markers import clean_mask, extract_bands
from .config import PANELS, CSV_EPOCHS


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


def extract_panel(rgb: np.ndarray, dark: np.ndarray, box: FluxBox) -> PanelResult:
    meta = PANELS[(box.figure, box.row, box.col)]
    calib = calibrate(dark, box)
    mm = clean_mask(rgb, box)
    upper, lower = extract_bands(mm, len(meta["bands"]))

    res = PanelResult(box.figure, box.row, box.col, meta["epoch"],
                      meta["offset"], meta["epoch"] in CSV_EPOCHS, calib)

    # bands[0] = I_C (optical, upper, circle); bands[1] = IR (lower, cross)
    bjd, flux = _pts_to_arrays(upper, calib, box)
    res.bands.append(BandData(meta["bands"][0], True, bjd, flux))
    if len(meta["bands"]) > 1:
        bjd, flux = _pts_to_arrays(lower, calib, box)
        res.bands.append(BandData(meta["bands"][1], False, bjd, flux))
    return res


def extract_figure(figure: str) -> list[PanelResult]:
    rgb = load_rgb(figure)
    dark = dark_mask(rgb)
    return [extract_panel(rgb, dark, b) for b in detect_flux_boxes(figure, rgb)]


def extract_all() -> list[PanelResult]:
    out = []
    for fig in ("f12", "f13"):
        out.extend(extract_figure(fig))
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
