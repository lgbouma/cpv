"""Detect data markers (filled-circle I_C and cross IR) inside a flux panel.

Strategy
--------
1. Build the black marker mask (model removed) restricted to the plot interior.
2. Exclude the legend (a light-gray framed box) and the epoch/axis text.
3. Estimate the marker radius from the mask.
4. Matched-filter the mask with a filled-disk kernel and an x-cross kernel;
   take local maxima with non-max suppression to get candidate centers.
5. Classify each candidate as circle vs cross from local patch features
   (fill fraction + disk/cross response), cross-checked by vertical band.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy import ndimage as ndi

from .image_io import gray
from .color_masks import marker_mask
from .panels import FluxBox


def _longest_run_ext(mask_1d: np.ndarray) -> tuple[int, int, int]:
    idx = np.flatnonzero(mask_1d)
    if idx.size == 0:
        return 0, 0, 0
    breaks = np.flatnonzero(np.diff(idx) > 1)
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks, [idx.size - 1]))
    lens = idx[ends] - idx[starts] + 1
    k = int(lens.argmax())
    return int(lens[k]), int(idx[starts[k]]), int(idx[ends[k]])


def detect_legend_box(rgb: np.ndarray, box: FluxBox, pad: int = 6,
                      min_run: int = 60):
    """Bounding box (local coords) of the legend frame, or None.

    The legend frame is a thin, neutral light-gray rounded rectangle; its top
    and bottom edges are the only long horizontal neutral-gray runs inside the
    plot. We locate those edges and take their pixel extent."""
    sub = rgb[box.y0:box.y1 + 1, box.x0:box.x1 + 1].astype(np.int16)
    g = gray(rgb[box.y0:box.y1 + 1, box.x0:box.x1 + 1])
    R, G, B = sub[..., 0], sub[..., 1], sub[..., 2]
    neutral = (np.abs(R - G) < 22) & (np.abs(R - B) < 22)
    framish = (g < 235) & (g > 110) & neutral
    H, W = g.shape
    edges = []
    for y in range(2, H - 2):
        L, a, c = _longest_run_ext(framish[y, 3:W - 3])
        if L > min_run:
            edges.append((y, a + 3, c + 3))
    if not edges:
        return None
    ys = [e[0] for e in edges]
    left = min(e[1] for e in edges)
    right = max(e[2] for e in edges)
    top, bot = min(ys), max(ys)
    return (max(0, left - pad), max(0, top - pad),
            min(W - 1, right + pad), min(H - 1, bot + pad))


def _disk(r: int) -> np.ndarray:
    yy, xx = np.ogrid[-r:r + 1, -r:r + 1]
    return (xx * xx + yy * yy) <= r * r


def detect_text_box(mask: np.ndarray, region=(0.58, 0.0, 1.0, 0.66),
                    min_w: int = 30, max_w: int = 185, max_h: int = 32,
                    pad: int = 4):
    """Bounding box of the epoch label (e.g. "2014 Dec 27"), a short text
    string at the bottom-left, below the data bands.

    After horizontal dilation the data bands merge into wide blobs (>~200 px)
    while the epoch label stays a narrow word (~80-140 px). We keep narrow,
    short, bottom-most blobs whose centroid sits in the lowest fifth of the
    panel (where no band lives), to avoid masking real data.
    """
    H, W = mask.shape
    y0 = int(region[0] * H); x0 = int(region[1] * W)
    y1 = int(region[2] * H); x1 = int(region[3] * W)
    sub = np.zeros_like(mask)
    sub[y0:y1, x0:x1] = mask[y0:y1, x0:x1]
    if sub.sum() < 25:
        return None
    d = ndi.binary_dilation(sub, structure=np.ones((3, 21)))  # merge chars
    lab, n = ndi.label(d)
    cands = []
    for i in range(1, n + 1):
        ys, xs = np.where(lab == i)
        w = xs.max() - xs.min(); h = ys.max() - ys.min()
        if min_w <= w <= max_w and h <= max_h and ys.mean() > 0.78 * H:
            cands.append((ys.mean(), xs.min(), ys.min(), xs.max(), ys.max()))
    if not cands:
        return None
    cands.sort(key=lambda c: -c[0])  # bottom-most first
    _, xa, ya, xb, yb = cands[0]
    return (max(0, xa - pad), max(0, ya - pad),
            min(W - 1, xb + pad), min(H - 1, yb + pad))


def clean_mask(rgb: np.ndarray, box: FluxBox, margin: int = 10) -> np.ndarray:
    """Marker mask with model, legend, epoch text, and frame margins removed.

    The margin clears the inward-pointing axis tick marks (~8 px) on all four
    spines; the data sits well inside this margin in every panel."""
    sub = rgb[box.y0:box.y1 + 1, box.x0:box.x1 + 1]
    mm = marker_mask(sub).copy()
    mm[:margin, :] = False; mm[-margin:, :] = False
    mm[:, :margin] = False; mm[:, -margin:] = False
    lb = detect_legend_box(rgb, box)
    if lb:
        mm[lb[1]:lb[3] + 1, lb[0]:lb[2] + 1] = False
    tb = detect_text_box(mm)
    if tb:
        mm[tb[1]:tb[3] + 1, tb[0]:tb[2] + 1] = False
    return mm


def _reject_outliers(pts, thresh_px: float = 22.0, win: int = 15):
    """Drop points far from a rolling-median trend (isolated strays from ticks
    or stray text); keeps gradual dip structure, which the median follows."""
    if len(pts) < 6:
        return pts
    arr = np.array(pts)
    order = np.argsort(arr[:, 0])
    arr = arr[order]
    ymed = ndi.median_filter(arr[:, 1], size=win, mode="nearest")
    keep = np.abs(arr[:, 1] - ymed) <= thresh_px
    return [tuple(p) for p in arr[keep]]


def _column_ys(mask: np.ndarray, x_lo: int, x_hi: int) -> np.ndarray:
    """Sorted y pixel positions of mask hits in columns [x_lo, x_hi)."""
    cols = mask[:, x_lo:x_hi]
    ys = np.where(cols.any(1))[0]
    return ys


def extract_bands(mask: np.ndarray, n_bands: int, xbin: int = 2,
                  min_pix: int = 2, gap_min: int = 14):
    """Extract per-column (x_px, y_px) points for each band.

    For two-band panels the optical (I_C, circles) curve is always the upper
    band and the IR (crosses) the lower; we split each column at its dominant
    y-gap. Columns with a single cluster (one band absent there) are assigned
    to the nearer band by continuity in a second pass.

    Returns (upper_pts, lower_pts) for n_bands==2, else (single_pts, []).
    """
    H, W = mask.shape
    xs_centers = list(range(0, W - xbin + 1, xbin))

    if n_bands == 1:
        pts = []
        for xc in xs_centers:
            ys = _column_ys(mask, xc, xc + xbin)
            if ys.size >= min_pix:
                pts.append((xc + (xbin - 1) / 2.0, float(ys.mean())))
        return _reject_outliers(pts), []

    clear_up, clear_lo, ambiguous = [], [], []
    for xc in xs_centers:
        ys = _column_ys(mask, xc, xc + xbin)
        if ys.size < min_pix:
            continue
        xcen = xc + (xbin - 1) / 2.0
        if ys.size == 1:
            ambiguous.append((xcen, float(ys[0])))
            continue
        gaps = np.diff(ys)
        k = int(gaps.argmax())
        if gaps[k] >= gap_min:
            clear_up.append((xcen, float(ys[:k + 1].mean())))
            clear_lo.append((xcen, float(ys[k + 1:].mean())))
        else:
            ambiguous.append((xcen, float(ys.mean())))

    # build trends from the clear splits to assign ambiguous columns
    up = np.array(clear_up) if clear_up else np.empty((0, 2))
    lo = np.array(clear_lo) if clear_lo else np.empty((0, 2))

    def _interp(trend, x):
        if trend.shape[0] == 0:
            return None
        return float(np.interp(x, trend[:, 0], trend[:, 1]))

    for xcen, ymean in ambiguous:
        yu = _interp(up, xcen)
        yl = _interp(lo, xcen)
        if yu is None and yl is None:
            continue
        if yu is None:
            clear_lo.append((xcen, ymean)); continue
        if yl is None:
            clear_up.append((xcen, ymean)); continue
        if abs(ymean - yu) <= abs(ymean - yl):
            clear_up.append((xcen, ymean))
        else:
            clear_lo.append((xcen, ymean))

    clear_up.sort(); clear_lo.sort()
    return _reject_outliers(clear_up), _reject_outliers(clear_lo)
