"""Reproduce Tanimoto-style panels from digitized data, plus QA overlays."""
from __future__ import annotations

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .image_io import load_rgb
from .panels import detect_flux_boxes
from .config import PANELS
from .extract import PanelResult

# panel grid shape per figure (rows, cols)
GRID = {"f12": (4, 3), "f13": (3, 3)}


def _panel_axis(results, figure, row, col):
    for r in results:
        if r.figure == figure and r.row == row and r.col == col:
            return r
    return None


def plot_comparison(figure: str, results: list[PanelResult], outpath: str):
    """Reproduce the figure's panel grid using only digitized points."""
    nrow, ncol = GRID[figure]
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.0 * nrow))
    axs = np.atleast_2d(axs)
    used = set()
    for r in results:
        if r.figure != figure:
            continue
        used.add((r.row, r.col))
        ax = axs[r.row, r.col]
        for bd in r.bands:
            if bd.bjd.size == 0:
                continue
            if bd.is_optical:
                ax.plot(bd.bjd, bd.flux, "o", ms=2.6, color="k",
                        label=_band_label(bd.band))
            else:
                ax.plot(bd.bjd, bd.flux, "x", ms=3.4, mew=0.7, color="k",
                        label=_band_label(bd.band))
        ax.text(0.03, 0.06, r.epoch, transform=ax.transAxes, fontsize=9)
        ax.set_xlabel(f"BJD $-$ {r.offset}", fontsize=9)
        ax.set_ylabel("Relative flux", fontsize=9)
        ax.set_xlim(*r.xlim)
        ax.set_ylim(*r.ylim)
        ax.legend(loc="lower right", fontsize=7, framealpha=0.9)
        ax.tick_params(labelsize=8)
    # hide unused axes
    for rr in range(nrow):
        for cc in range(ncol):
            if (rr, cc) not in used:
                axs[rr, cc].axis("off")
    fig.suptitle(f"Digitized reproduction of {figure} (PTFO 8-8695, "
                 f"Tanimoto et al. 2020)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(outpath, dpi=130)
    plt.close(fig)


def _band_label(band: str) -> str:
    mapping = {"K_s": "$K_s$", "J": "$J$", "H": "$H$", "I_C": "$I_C$"}
    return mapping.get(band, band) + " band"


def plot_side_by_side(results: list[PanelResult], outpath: str,
                      only_csv: bool = True):
    """Original panel crop (left) vs digitized reproduction (right), one row
    per epoch. Defaults to the epochs overlapping Table 3."""
    sel = [r for r in results if (r.in_csv or not only_csv)]
    sel.sort(key=lambda r: (r.figure, r.row, r.col))
    n = len(sel)
    fig, axs = plt.subplots(n, 2, figsize=(12, 2.7 * n))
    axs = np.atleast_2d(axs)
    rgb_cache = {}
    boxes_cache = {}
    for i, r in enumerate(sel):
        if r.figure not in rgb_cache:
            rgb_cache[r.figure] = load_rgb(r.figure)
            boxes_cache[r.figure] = {(b.row, b.col): b
                                     for b in detect_flux_boxes(r.figure, rgb_cache[r.figure])}
        b = boxes_cache[r.figure][(r.row, r.col)]
        sub = rgb_cache[r.figure][b.y0:b.y1 + 1, b.x0:b.x1 + 1]
        # Draw the original crop in data coordinates (extent + aspect="auto")
        # so its axes span exactly the panel's xlim/ylim. Both columns then
        # share the same x-width and the dips line up horizontally, making the
        # digitization quality easy to read off by eye.
        axl = axs[i, 0]
        axl.imshow(sub, extent=[r.xlim[0], r.xlim[1], r.ylim[0], r.ylim[1]],
                   aspect="auto")
        axl.set_xlim(*r.xlim)
        axl.set_ylim(*r.ylim)
        axl.set_title(f"original: {r.epoch}", fontsize=9)
        axl.set_xlabel(f"BJD $-$ {r.offset}", fontsize=8)
        axl.set_ylabel("Relative flux", fontsize=8)
        axl.tick_params(labelsize=7)
        axr = axs[i, 1]
        for bd in r.bands:
            if bd.bjd.size == 0:
                continue
            mk = "o" if bd.is_optical else "x"
            axr.plot(bd.bjd, bd.flux, mk, ms=2.8, mew=0.7, color="k",
                     label=_band_label(bd.band))
        axr.set_title(f"digitized: {r.epoch}", fontsize=9)
        axr.set_xlabel(f"BJD $-$ {r.offset}", fontsize=8)
        axr.set_ylabel("Relative flux", fontsize=8)
        axr.set_xlim(*r.xlim)
        axr.set_ylim(*r.ylim)
        axr.legend(loc="lower right", fontsize=7)
        axr.tick_params(labelsize=7)
    fig.suptitle("PTFO 8-8695 (Tanimoto et al. 2020): original vs digitized",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.995))
    fig.savefig(outpath, dpi=115)
    plt.close(fig)


def _band_pixels(bd, calib, box):
    """Invert the calibration to get a band's box-local (x_px, y_px), so an
    overlay shows exactly the points that were extracted (any method)."""
    xpx = (bd.bjd - calib.x_intercept) / calib.x_slope - box.x0
    ypx = (bd.flux - calib.y_intercept) / calib.y_slope - box.y0
    return xpx, ypx


def plot_diagnostic_overlays(figure: str, results: list[PanelResult],
                             outdir: str):
    """Combined sheet of the extracted points drawn on the original crops, for
    visual QA. Points come from `results`, so the overlay reflects whichever
    extraction method produced them."""
    rgb = load_rgb(figure)
    boxes = {(b.row, b.col): b for b in detect_flux_boxes(figure, rgb)}
    nrow, ncol = GRID[figure]
    sheet, sheet_axs = plt.subplots(nrow, ncol, figsize=(5.0 * ncol, 3.0 * nrow))
    sheet_axs = np.atleast_2d(sheet_axs)
    used = set()
    for r in results:
        if r.figure != figure:
            continue
        b = boxes[(r.row, r.col)]
        sub = rgb[b.y0:b.y1 + 1, b.x0:b.x1 + 1]
        ax = sheet_axs[r.row, r.col]
        ax.imshow(sub)
        for bd in r.bands:
            if bd.bjd.size == 0:
                continue
            xpx, ypx = _band_pixels(bd, r.calib, b)
            if bd.is_optical:
                ax.scatter(xpx, ypx, s=5, c="blue", marker="o", lw=0)
            else:
                ax.scatter(xpx, ypx, s=13, c="lime", marker="x", lw=0.7)
        ax.set_title(f"{r.epoch}  ({'+'.join(bd.band for bd in r.bands)})",
                     fontsize=8)
        ax.axis("off")
        used.add((r.row, r.col))
    for rr in range(nrow):
        for cc in range(ncol):
            if (rr, cc) not in used:
                sheet_axs[rr, cc].axis("off")
    sheet.suptitle(f"{figure}: digitized points (blue=I_C circles, "
                   f"green=IR crosses) on originals", fontsize=12)
    sheet.tight_layout(rect=(0, 0, 1, 0.98))
    sheet.savefig(os.path.join(outdir, f"overlay_{figure}.png"), dpi=110)
    plt.close(sheet)


def plot_method_comparison(results_by_method: dict, outpath: str,
                           only_csv: bool = True):
    """One row per epoch: original crop vs each method's digitized panel, all
    sharing the panel's true xlim/ylim. The key artifact for judging how the
    two extraction methods differ."""
    names = list(results_by_method)
    base = [r for r in results_by_method[names[0]] if (r.in_csv or not only_csv)]
    base.sort(key=lambda r: (r.figure, r.row, r.col))
    n = len(base)
    ncol = 1 + len(names)
    fig, axs = plt.subplots(n, ncol, figsize=(3.9 * ncol, 2.5 * n))
    axs = np.atleast_2d(axs)
    rgb_cache, boxes_cache = {}, {}
    for i, r0 in enumerate(base):
        if r0.figure not in rgb_cache:
            rgb_cache[r0.figure] = load_rgb(r0.figure)
            boxes_cache[r0.figure] = {(b.row, b.col): b for b in
                                      detect_flux_boxes(r0.figure, rgb_cache[r0.figure])}
        b = boxes_cache[r0.figure][(r0.row, r0.col)]
        sub = rgb_cache[r0.figure][b.y0:b.y1 + 1, b.x0:b.x1 + 1]
        axl = axs[i, 0]
        axl.imshow(sub, extent=[r0.xlim[0], r0.xlim[1], r0.ylim[0], r0.ylim[1]],
                   aspect="auto")
        axl.set_xlim(*r0.xlim); axl.set_ylim(*r0.ylim)
        axl.set_title(f"original: {r0.epoch}", fontsize=8)
        axl.set_ylabel("Rel. flux", fontsize=7)
        axl.tick_params(labelsize=6)
        for j, name in enumerate(names):
            r = _panel_axis(results_by_method[name], r0.figure, r0.row, r0.col)
            ax = axs[i, 1 + j]
            npts = 0
            for bd in r.bands:
                if bd.bjd.size == 0:
                    continue
                npts += bd.bjd.size
                mk = "o" if bd.is_optical else "x"
                ax.plot(bd.bjd, bd.flux, mk, ms=2.6, mew=0.6, color="k")
            ax.set_xlim(*r0.xlim); ax.set_ylim(*r0.ylim)
            ax.set_title(f"{name}: {npts} pts", fontsize=8)
            ax.tick_params(labelsize=6)
    fig.suptitle("PTFO 8-8695 (Tanimoto+2020): original vs extraction methods",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.995))
    fig.savefig(outpath, dpi=120)
    plt.close(fig)
