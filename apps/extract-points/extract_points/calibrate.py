"""Pixel <-> data calibration from detected major tick marks.

Major ticks point inward from the spines. We read the x-tick pixel columns
just inside the flux box's bottom (BJD) axis and the y-tick pixel rows just
inside the left (relative flux) axis, then fit a linear pixel->data transform.
The BJD tick values come from config (`x_first` + 0.1*i); the flux tick values
are 1.0, 0.9, ... (top = 1.0).
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .config import PANELS, X_STEP, Y_TOP, Y_STEP
from .panels import FluxBox


@dataclass
class Calibration:
    # data = slope * pixel + intercept
    x_slope: float
    x_intercept: float
    y_slope: float
    y_intercept: float
    x_ticks_px: list[float]
    y_ticks_px: list[float]
    x_r2: float
    y_r2: float

    def px_to_bjd(self, px) -> np.ndarray:
        return self.x_slope * np.asarray(px, float) + self.x_intercept

    def px_to_flux(self, py) -> np.ndarray:
        return self.y_slope * np.asarray(py, float) + self.y_intercept


def _longest_run(mask_1d: np.ndarray) -> int:
    idx = np.flatnonzero(mask_1d)
    if idx.size == 0:
        return 0
    breaks = np.flatnonzero(np.diff(idx) > 1)
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks, [idx.size - 1]))
    return int((idx[ends] - idx[starts] + 1).max())


def _cluster_centers(flags: np.ndarray, base: int) -> list[float]:
    """Centers (in absolute pixel coords) of consecutive True runs in `flags`."""
    idx = np.flatnonzero(flags)
    if idx.size == 0:
        return []
    breaks = np.flatnonzero(np.diff(idx) > 1)
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks, [idx.size - 1]))
    return [base + 0.5 * (idx[s] + idx[e]) for s, e in zip(starts, ends)]


def find_bjd_axis_y(dark: np.ndarray, box: FluxBox, search: int = 260) -> int:
    """Pixel row of the labeled BJD axis: the bottom spine of the residual
    subplot below the flux box. It is the lowest full-width horizontal line
    before the inter-row whitespace. The shared x-transform matches the flux
    subplot exactly, but ticks are reliably present only on this axis."""
    width = box.x1 - box.x0
    y_end = min(dark.shape[0] - 1, box.y1 + search)
    hits = []
    for y in range(box.y1 + 5, y_end):
        run = _longest_run(dark[y, box.x0:box.x1 + 1])
        if run >= 0.85 * width:
            hits.append(y)
    if not hits:
        return box.y1  # fallback: flux-box bottom
    # cluster consecutive rows into lines, take the bottom-most line center
    hits = np.array(hits)
    breaks = np.flatnonzero(np.diff(hits) > 4)
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks, [hits.size - 1]))
    centers = [int(round(0.5 * (hits[s] + hits[e]))) for s, e in zip(starts, ends)]
    return max(centers)


def _inward_run_lengths(dark: np.ndarray, axis_y: int, x0: int, x1: int,
                        maxrun: int = 16) -> np.ndarray:
    """For each column, the number of consecutive dark pixels going upward
    (inward) from just above the axis. Tick marks show as nonzero runs; major
    ticks are longer (~4-5 px) than minor ticks (<=3 px)."""
    cols = np.arange(x0, x1 + 1)
    L = np.zeros(cols.size, int)
    for j, x in enumerate(cols):
        n = 0
        for y in range(axis_y - 2, axis_y - 2 - maxrun, -1):
            if y < 0 or not dark[y, x]:
                break
            n += 1
        L[j] = n
    return L


def detect_x_ticks(dark: np.ndarray, box: FluxBox) -> list[float]:
    """Major tick columns on the BJD axis (residual subplot bottom spine).

    Major ticks point inward (up). They are distinguished from the shorter
    minor ticks by their inward run length: we cluster dark columns into ticks,
    measure each tick's length, and keep the long (major) ones."""
    axis_y = find_bjd_axis_y(dark, box)
    L = _inward_run_lengths(dark, axis_y, box.x0, box.x1)
    cand = L >= 2  # any tick (major or minor)
    # cluster adjacent candidate columns into ticks; record center and length
    idx = np.flatnonzero(cand)
    if idx.size == 0:
        return []
    breaks = np.flatnonzero(np.diff(idx) > 1)
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks, [idx.size - 1]))
    ticks = []  # (center_px, length)
    for s, e in zip(starts, ends):
        seg = idx[s:e + 1]
        center = box.x0 + 0.5 * (seg[0] + seg[-1])
        ticks.append((center, int(L[seg].max())))
    # drop the spine corners (very tall runs at x0/x1)
    ticks = [(c, l) for c, l in ticks if abs(c - box.x0) > 3 and abs(c - box.x1) > 3]
    if not ticks:
        return []
    max_len = max(l for _, l in ticks)
    thr = max(4, max_len - 1)  # majors are the longest; minors are <=3 px
    majors = sorted(c for c, l in ticks if l >= thr)
    if len(majors) >= 2:
        s = np.median(np.diff(majors))  # major spacing
        # A major tick may coincide with a spine (xlim starts/ends on a tick),
        # where it merges with the vertical spine line and is missed. Recover
        # it by extrapolation if the predicted position lands on a spine.
        if abs((majors[0] - s) - box.x0) <= 4:
            majors.insert(0, majors[0] - s)
        if abs((majors[-1] + s) - box.x1) <= 4:
            majors.append(majors[-1] + s)
    return sorted(majors)


def detect_y_ticks(dark: np.ndarray, box: FluxBox, tlen: int = 8,
                   frac: float = 0.55) -> list[float]:
    """Major tick rows just inside the left (relative flux) axis.

    Major ticks point inward (right) from the left spine."""
    x_lo = box.x0 + 1
    x_hi = x_lo + tlen
    band = dark[box.y0:box.y1 + 1, x_lo:x_hi]
    rowcount = band.sum(1)
    flags = rowcount >= frac * (x_hi - x_lo)
    rows = _cluster_centers(flags, box.y0)
    rows = [r for r in rows if abs(r - box.y0) > 3 and abs(r - box.y1) > 3]
    return sorted(rows)


def _linfit(px: np.ndarray, val: np.ndarray) -> tuple[float, float, float]:
    """Least-squares val = slope*px + intercept; returns (slope, intercept, R2)."""
    px = np.asarray(px, float)
    val = np.asarray(val, float)
    A = np.vstack([px, np.ones_like(px)]).T
    (slope, intercept), *_ = np.linalg.lstsq(A, val, rcond=None)
    pred = slope * px + intercept
    ss_res = float(np.sum((val - pred) ** 2))
    ss_tot = float(np.sum((val - val.mean()) ** 2)) or 1.0
    r2 = 1.0 - ss_res / ss_tot
    return float(slope), float(intercept), r2


def _largest_uniform_subset(px: list[float], step_tol: float = 0.18) -> list[float]:
    """Keep ticks consistent with a single uniform spacing (robust to a stray
    detection). Uses the median spacing of sorted ticks and drops outliers."""
    if len(px) <= 2:
        return px
    px = sorted(px)
    diffs = np.diff(px)
    med = np.median(diffs)
    keep = [px[0]]
    for a, b in zip(px[:-1], px[1:]):
        # accept if gap is ~ an integer multiple of the median spacing
        ratio = (b - a) / med
        if abs(ratio - round(ratio)) <= step_tol and round(ratio) >= 1:
            keep.append(b)
    return keep


def calibrate(dark: np.ndarray, box: FluxBox) -> Calibration:
    meta = PANELS[(box.figure, box.row, box.col)]

    xpx = _largest_uniform_subset(detect_x_ticks(dark, box))
    xpx = sorted(xpx)
    xvals = meta["x_first"] + X_STEP * np.arange(len(xpx))
    xs, xi, xr2 = _linfit(np.array(xpx), xvals)

    ypx = sorted(detect_y_ticks(dark, box))  # top -> bottom (ascending pixel)
    yvals = Y_TOP - Y_STEP * np.arange(len(ypx))  # 1.0, 0.9, ...
    ys, yi, yr2 = _linfit(np.array(ypx), yvals)

    return Calibration(xs, xi, ys, yi, xpx, ypx, xr2, yr2)
