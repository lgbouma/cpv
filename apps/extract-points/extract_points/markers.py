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
from .color_masks import marker_mask, black_mask, red_mask
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


def reclaim_occluded_red(sub: np.ndarray, reach: int = 14) -> np.ndarray:
    """Black markers with red-model-occluded interiors filled back in.

    Tanimoto overlay a red model curve on the black markers; ``marker_mask``
    removes red, which punches holes through marker centers wherever the model
    crosses them (worst inside the dip, where the model lives). We reclaim a red
    pixel iff it is vertically sandwiched between black within ``reach`` px
    (~one marker diameter): such pixels are occluded marker interiors. The bare
    model curve in data gaps and between the two bands has no flanking black, so
    it is left out."""
    bk = black_mask(sub)
    rd = red_mask(sub)
    above = np.zeros_like(bk)
    below = np.zeros_like(bk)
    for k in range(1, reach + 1):
        above[k:, :] |= bk[:-k, :]   # black exists within reach px above
        below[:-k, :] |= bk[k:, :]   # black exists within reach px below
    return bk | (rd & above & below)


def clean_mask(rgb: np.ndarray, box: FluxBox, margin: int = 10,
               red_fill: bool = False, return_black: bool = False):
    """Marker mask with model, legend, epoch text, and frame margins removed.

    The margin clears the inward-pointing axis tick marks (~8 px) on all four
    spines; the data sits well inside this margin in every panel. With
    ``red_fill`` the red-occluded marker interiors are restored (see
    ``reclaim_occluded_red``) instead of being dropped.

    With ``return_black`` we also return the genuine-black subset of the cleaned
    mask (``mm & black_mask``). For ``red_fill`` this is the real marker pixels
    *without* the reclaimed model curve, so a reducer can detect markers from the
    filled mask but read each point's y (flux) off the black markers only --
    keeping the extraction from tracing Tanimoto's overlaid red model."""
    sub = rgb[box.y0:box.y1 + 1, box.x0:box.x1 + 1]
    mm = (reclaim_occluded_red(sub) if red_fill else marker_mask(sub)).copy()
    mm[:margin, :] = False; mm[-margin:, :] = False
    mm[:, :margin] = False; mm[:, -margin:] = False
    lb = detect_legend_box(rgb, box)
    if lb:
        mm[lb[1]:lb[3] + 1, lb[0]:lb[2] + 1] = False
    tb = detect_text_box(mm)
    if tb:
        mm[tb[1]:tb[3] + 1, tb[0]:tb[2] + 1] = False
    if return_black:
        return mm, mm & black_mask(sub)
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


# ---------------------------------------------------------------------------
# Band separation into per-band PIXEL masks (shared by both extraction methods)
# ---------------------------------------------------------------------------
def split_band_masks(mask: np.ndarray, n_bands: int, xbin: int = 2,
                     min_pix: int = 2, gap_min: int = 14) -> list[np.ndarray]:
    """Assign every marker pixel to a band, returning one boolean mask per band.

    Reuses the per-column y-gap split that ``extract_bands`` uses to separate
    the upper (optical I_C) band from the lower (IR) band, but instead of
    collapsing each column to a single mean point it keeps the individual
    pixels. This lets a downstream reducer (matched-filter marker detection or
    cadence binning) work from the real marker pixels.

    Returns ``[band0]`` for 1-band panels, else ``[upper, lower]`` (band0 is
    always the optical/circle band).
    """
    H, W = mask.shape
    if n_bands == 1:
        return [mask.copy()]

    band0 = np.zeros_like(mask)
    band1 = np.zeros_like(mask)
    ygrid = np.arange(H)
    xs_centers = list(range(0, W - xbin + 1, xbin))

    cols = []  # (xc, xcen, ysplit, mean) deferred to the assignment pass
    clear_up, clear_lo = [], []
    for xc in xs_centers:
        ys = _column_ys(mask, xc, xc + xbin)
        if ys.size < min_pix:
            continue
        xcen = xc + (xbin - 1) / 2.0
        if ys.size == 1:
            cols.append((xc, xcen, None, float(ys[0])))
            continue
        gaps = np.diff(ys)
        k = int(gaps.argmax())
        if gaps[k] >= gap_min:
            ysplit = 0.5 * (ys[k] + ys[k + 1])
            cols.append((xc, xcen, ysplit, None))
            clear_up.append((xcen, float(ys[:k + 1].mean())))
            clear_lo.append((xcen, float(ys[k + 1:].mean())))
        else:
            cols.append((xc, xcen, None, float(ys.mean())))

    up = np.array(clear_up) if clear_up else np.empty((0, 2))
    lo = np.array(clear_lo) if clear_lo else np.empty((0, 2))

    def _interp(trend, x):
        if trend.shape[0] == 0:
            return None
        return float(np.interp(x, trend[:, 0], trend[:, 1]))

    for xc, xcen, ysplit, mean in cols:
        seg = slice(xc, xc + xbin)
        sub = mask[:, seg]
        if ysplit is not None:
            up_sel = (ygrid <= ysplit)[:, None]
            band0[:, seg] = sub & up_sel
            band1[:, seg] = sub & ~up_sel
        else:
            yu = _interp(up, xcen)
            yl = _interp(lo, xcen)
            if yu is None and yl is None:
                continue
            to_upper = (yl is None) or (yu is not None and
                                        abs(mean - yu) <= abs(mean - yl))
            if to_upper:
                band0[:, seg] |= sub
            else:
                band1[:, seg] |= sub
    return [band0, band1]


def _estimate_geom(band_mask: np.ndarray) -> tuple[int, int]:
    """Rough (radius, pitch) in px from the band's per-column vertical extent.

    The median per-column extent reflects ~2 overlapping markers, so a marker
    radius ~ median/3 and a marker pitch ~ one diameter is a reasonable, if
    approximate, estimate. Both are clamped to sane ranges and used only to set
    the matched-filter kernel size and the cadence bin width."""
    H, W = band_mask.shape
    ext = []
    for x in range(W):
        ys = np.flatnonzero(band_mask[:, x])
        if ys.size:
            ext.append(ys.max() - ys.min() + 1)
    if not ext:
        return 4, 8
    med = float(np.median(ext))
    r = int(np.clip(round(med / 3.0), 3, 8))
    # The ribbon runs ~2 markers thick (median extent ~2 diameters), i.e. the
    # true marker pitch is ~one radius, not one diameter; sampling at the radius
    # recovers ~one point per real marker without the per-column oversampling.
    pitch = max(3, r)
    return r, pitch


def _comb_centers(cand: np.ndarray, det_mask: np.ndarray,
                  pos_mask: np.ndarray,
                  pitch: int) -> list[tuple[float, float]]:
    """Turn thresholded matched-filter cores into one center per ~marker.

    Each connected core region is split into ``round(width / pitch)`` evenly
    spaced x-positions (so a solid, fully-overlapping ribbon yields a comb at
    the marker pitch while a resolved/isolated marker yields a single point).

    The x position and the comb geometry come from ``det_mask`` (the detection
    mask); the y (flux) is the median of ``pos_mask`` pixels in the same
    half-pitch window. For per_marker/cadence the two masks are identical. For
    red_aware ``det_mask`` is the red-filled mask (so markers detect cleanly
    even where the model punches holes) while ``pos_mask`` is the genuine-black
    subset -- so the flux is read off the real markers, not the reclaimed model
    curve. If a window has no black pixel we fall back to the detection mask."""
    lab, n = ndi.label(cand)
    hw = max(2, pitch // 2)
    centers = []
    for i in range(1, n + 1):
        xs = np.where(lab == i)[1]
        xa, xb = int(xs.min()), int(xs.max())
        k = max(1, int(round((xb - xa + 1) / float(pitch))))
        xpos = [0.5 * (xa + xb)] if k == 1 else list(np.linspace(xa, xb, k))
        for xp in xpos:
            xi = int(round(xp))
            lo_x = max(0, xi - hw)
            dy, dx = np.nonzero(det_mask[:, lo_x:xi + hw + 1])
            if dx.size == 0:
                continue
            py, _ = np.nonzero(pos_mask[:, lo_x:xi + hw + 1])
            yvals = py if py.size else dy
            centers.append((float(np.mean(dx)) + lo_x, float(np.median(yvals))))
    centers.sort()
    return centers


def detect_markers_matched(band_mask: np.ndarray, pos_mask: np.ndarray = None,
                           thr_frac: float = 0.45) -> list[tuple[float, float]]:
    """Option 1: per-marker detection by matched filtering with a disk kernel.

    Correlate the band mask with a filled-disk kernel (sized from the data),
    keep local maxima of the normalized response above ``thr_frac`` as marker
    cores, then resolve one centroid per marker via ``_comb_centers``. The disk
    kernel localizes both circle (I_C) and cross (IR) clusters; the actual
    center is the centroid of the real pixels, so it is shape-agnostic.

    Detection and the marker-pitch comb run on ``band_mask``; the y (flux) of
    each point is read from ``pos_mask`` if given (else ``band_mask``). red_aware
    passes the genuine-black mask here so it detects from the red-filled mask but
    measures flux off the black markers -- see ``_comb_centers``.

    Returns box-local (x_px, y_px) centers."""
    if band_mask.sum() == 0:
        return []
    if pos_mask is None:
        pos_mask = band_mask
    r, pitch = _estimate_geom(band_mask)
    kern = _disk(r).astype(float)
    resp = ndi.correlate(band_mask.astype(float), kern, mode="constant")
    if resp.max() <= 0:
        return []
    rn = resp / resp.max()
    # Core = where a full marker disk fits (high response). Each connected core
    # is one resolved marker (isolated) or a merged run; _comb_centers emits a
    # marker-pitch comb across it, placing each point at a local pixel centroid.
    cand = rn >= thr_frac
    return _reject_outliers(_comb_centers(cand, band_mask, pos_mask, pitch))


def bin_cadence(band_mask: np.ndarray,
                pos_mask: np.ndarray = None) -> list[tuple[float, float]]:
    """Option 2: cadence-matched binning.

    Bin x at the estimated marker pitch and take the MEDIAN y of the band
    pixels in each bin (x = mean of their columns). Robust and simple; collapses
    the per-column oversampling and the inter-marker centroid wander.

    ``pos_mask`` (the genuine-black subset, if given) supplies the y values while
    ``band_mask`` gates the bin and sets x -- mirroring ``detect_markers_matched``
    so flux is read off the real markers. For per_marker/cadence the two coincide.

    Returns box-local (x_px, y_px) points."""
    if band_mask.sum() == 0:
        return []
    if pos_mask is None:
        pos_mask = band_mask
    H, W = band_mask.shape
    _, pitch = _estimate_geom(band_mask)
    pts = []
    for xa in range(0, W, pitch):
        dy, dx = np.nonzero(band_mask[:, xa:xa + pitch])
        if dy.size < 2:
            continue
        py, _ = np.nonzero(pos_mask[:, xa:xa + pitch])
        yvals = py if py.size else dy
        pts.append((float(np.mean(dx)) + xa, float(np.median(yvals))))
    pts.sort()
    return _reject_outliers(pts)
