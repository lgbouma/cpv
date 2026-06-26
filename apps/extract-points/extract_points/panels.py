"""Locate axes frames and isolate the flux subpanel of each panel.

The figures are regular grids. Each panel is a tall 'Relative flux' box
stacked over a short 'Residuals' box that share left/right spines. We detect
spine lines as long continuous runs of dark pixels, reconstruct the grid, and
return the flux box (the tall one) for every panel.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .image_io import dark_mask


@dataclass
class FluxBox:
    figure: str
    row: int
    col: int
    x0: int  # left spine pixel (inclusive)
    x1: int  # right spine pixel
    y0: int  # top spine pixel
    y1: int  # bottom spine pixel (flux/residual boundary, top of residuals)

    @property
    def index(self) -> int:
        return self._index

    def __post_init__(self):
        self._index = -1


def _longest_run_len(mask_1d: np.ndarray) -> int:
    """Length of the longest run of True in a 1-D boolean array."""
    if not mask_1d.any():
        return 0
    # positions where value changes
    idx = np.flatnonzero(mask_1d)
    if idx.size == 0:
        return 0
    # split into consecutive groups
    breaks = np.flatnonzero(np.diff(idx) > 1)
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks, [idx.size - 1]))
    runs = idx[ends] - idx[starts] + 1
    return int(runs.max())


def _line_positions(dark: np.ndarray, axis: int, min_len: int) -> list[tuple[int, int]]:
    """Find spine lines along `axis`.

    axis=0 -> horizontal lines (scan each row); returns clusters of row indices.
    axis=1 -> vertical lines (scan each column); returns clusters of col indices.
    Returns list of (center, thickness) pixel positions of detected lines.
    """
    n = dark.shape[axis]
    runs = np.empty(n, dtype=int)
    for i in range(n):
        line = dark[i, :] if axis == 0 else dark[:, i]
        runs[i] = _longest_run_len(line)
    is_line = runs >= min_len
    # cluster consecutive True indices into single lines
    out = []
    idx = np.flatnonzero(is_line)
    if idx.size == 0:
        return out
    breaks = np.flatnonzero(np.diff(idx) > 3)
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks, [idx.size - 1]))
    for s, e in zip(starts, ends):
        seg = idx[s:e + 1]
        # weight center by run length
        w = runs[seg].astype(float)
        center = int(round((seg * w).sum() / w.sum()))
        out.append((center, int(seg[-1] - seg[0] + 1)))
    return out


def _cluster_1d(values: list[int], gap: int) -> list[list[int]]:
    """Group sorted ints into clusters split where consecutive gap > `gap`."""
    if not values:
        return []
    vals = sorted(values)
    groups = [[vals[0]]]
    for v in vals[1:]:
        if v - groups[-1][-1] > gap:
            groups.append([v])
        else:
            groups[-1].append(v)
    return groups


def detect_flux_boxes(figure: str, rgb: np.ndarray) -> list[FluxBox]:
    """Detect every flux subpanel box in a figure."""
    dark = dark_mask(rgb)
    H, W = dark.shape

    # Vertical spines: long vertical dark runs (>= ~1/8 image height).
    vlines = _line_positions(dark, axis=1, min_len=H // 12)
    vx = [c for c, _ in vlines]
    # Horizontal spines: long horizontal runs (>= ~1/6 image width).
    hlines = _line_positions(dark, axis=0, min_len=W // 6)
    hy = [c for c, _ in hlines]

    # Column groups: vertical spines come in (left, right) pairs per column.
    # Cluster x positions that are close (within ~1.5% of width) into single
    # spine, then pair consecutive spines into columns.
    vclust = _cluster_1d(vx, gap=max(4, W // 200))
    vpos = [int(round(np.mean(g))) for g in vclust]
    vpos.sort()
    # pair: spines come as L,R,L,R,... -> consecutive pairs whose width is a
    # large fraction of inter-column spacing.
    columns = []
    i = 0
    while i + 1 < len(vpos):
        x0, x1 = vpos[i], vpos[i + 1]
        if x1 - x0 > W // 12:  # plausible panel width
            columns.append((x0, x1))
            i += 2
        else:
            i += 1
    columns.sort()

    # Horizontal spines: cluster nearby rows into single lines, then the flux
    # box of each panel-row is the consecutive spine pair whose vertical gap is
    # large (the tall 'Relative flux' subplot). Residual boxes (~55-75 px) and
    # inter-row whitespace (~170 px) are smaller; flux boxes are ~265-275 px.
    hclust = _cluster_1d(hy, gap=3)
    hpos = sorted(int(round(np.mean(g))) for g in hclust)
    gaps = [b - a for a, b in zip(hpos[:-1], hpos[1:])]
    max_gap = max(gaps) if gaps else 0

    flux_rows = []  # (top, bottom)
    for a, b in zip(hpos[:-1], hpos[1:]):
        gap = b - a
        if gap > 0.72 * max_gap and gap > 120:  # tall flux subplot
            flux_rows.append((a, b))
    flux_rows.sort()

    boxes: list[FluxBox] = []
    for r, (flux_top, flux_bot) in enumerate(flux_rows):
        for c, (x0, x1) in enumerate(columns):
            # require the flux box actually has marker content (some panels in
            # the last row are absent, e.g. f13 row3 col2)
            sub = dark[flux_top + 3:flux_bot - 3, x0 + 3:x1 - 3]
            if sub.size == 0 or sub.mean() < 0.01:
                continue
            boxes.append(FluxBox(figure, r, c, x0, x1, flux_top, flux_bot))

    boxes.sort(key=lambda b: (b.row, b.col))
    for k, b in enumerate(boxes):
        b._index = k
    return boxes
