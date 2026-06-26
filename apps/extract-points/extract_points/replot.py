"""Reproduce Tanimoto-style panels from digitized data, plus QA overlays."""
from __future__ import annotations

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .image_io import load_rgb
from .panels import detect_flux_boxes
from .markers import clean_mask, extract_bands
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
        axs[i, 0].imshow(sub)
        axs[i, 0].set_title(f"original: {r.epoch}", fontsize=9)
        axs[i, 0].axis("off")
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
        axr.legend(loc="lower right", fontsize=7)
        axr.tick_params(labelsize=7)
    fig.suptitle("PTFO 8-8695 (Tanimoto et al. 2020): original vs digitized",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.995))
    fig.savefig(outpath, dpi=115)
    plt.close(fig)


def plot_diagnostic_overlays(figure: str, outdir: str):
    """Save one overlay per panel (detected points on the original crop) and a
    combined sheet, for visual QA of the extraction."""
    rgb = load_rgb(figure)
    boxes = detect_flux_boxes(figure, rgb)
    nrow, ncol = GRID[figure]
    sheet, sheet_axs = plt.subplots(nrow, ncol, figsize=(5.0 * ncol, 3.0 * nrow))
    sheet_axs = np.atleast_2d(sheet_axs)
    used = set()
    for b in boxes:
        meta = PANELS[(figure, b.row, b.col)]
        sub = rgb[b.y0:b.y1 + 1, b.x0:b.x1 + 1]
        mm = clean_mask(rgb, b)
        up, lo = extract_bands(mm, len(meta["bands"]))
        for ax in (sheet_axs[b.row, b.col],):
            ax.imshow(sub)
            if up:
                a = np.array(up)
                ax.scatter(a[:, 0], a[:, 1], s=5, c="blue", marker="o", lw=0)
            if lo:
                a = np.array(lo)
                ax.scatter(a[:, 0], a[:, 1], s=13, c="lime", marker="x", lw=0.7)
            ax.set_title(f"{meta['epoch']}  ({'+'.join(meta['bands'])})",
                         fontsize=8)
            ax.axis("off")
        used.add((b.row, b.col))
    for rr in range(nrow):
        for cc in range(ncol):
            if (rr, cc) not in used:
                sheet_axs[rr, cc].axis("off")
    sheet.suptitle(f"{figure}: digitized points (blue=I_C circles, "
                   f"green=IR crosses) on originals", fontsize=12)
    sheet.tight_layout(rect=(0, 0, 1, 0.98))
    sheet.savefig(os.path.join(outdir, f"overlay_{figure}.png"), dpi=110)
    plt.close(sheet)
