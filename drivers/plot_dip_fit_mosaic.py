"""
Publication mosaic of every fitted dip epoch that enters the extinction-law
(beta) fit table (papers/Bouma_2026_cgcd/tables/tab_fitresults_rows.tex, whose
depths come from tables/dipdepths_thiswork_and_literature.csv).

There are 20 per-epoch table rows; TIC 300651846's three same-night epochs are
combined into a single phased panel, giving 18 panels split across two 3x3
figures (9 panels each).

Each panel stacks the bands observed at that epoch (bluest on top). For every
band it reproduces the fit-dips overview: the preferred (trapezoid) model, the
baseline with no dip, the +/-1 sigma depth band (the point of the figure), and
the data points -- black for in-/out-of-dip, gray for the excluded ("flare")
points. Bands are colored following the conventions of
plot_lp12502_multicolor_phase / plot_tic262400835_multicolor_phase /
plot_fourstar_vs_tess: B blue, g green, r red, z purple, Tierras & i orange,
TESS black (time-binned) over a gray raw underplot, and reddening to
maroon/brown in the near-IR.

Model curves/depths are pulled from the fit-dips registry (the joint
shared-duration dips when present, else the preferred single-band dips), exactly
as fit_dips.export_csv builds the depth table -- so the panels match the table.

Run (from repo root):
    apps/fit-dips is added to sys.path automatically.
    python drivers/plot_dip_fit_mosaic.py
"""
import csv
import os
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from aesthetic.plot import set_style, savefig

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "apps", "fit-dips"))

from fit_dips import registry, loaders, export_csv          # noqa: E402
from fit_dips.models import baseline_value, dip_model, parse_model_name  # noqa: E402

OUTDIR = os.path.join(REPO, "results", "dip_mosaic")
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Panel definition: TeX-table order, with the display title requested for each
# panel. Bands for each epoch are pulled from the depths CSV (grouped by
# star+date_ut), so the mosaic tracks whatever is in the beta-fit table.
# ---------------------------------------------------------------------------
PANELS = [
    ("TIC 262400835", "2020 Dec 13", "TIC 262400835, 2020 Dec 13"),
    ("TIC 262400835", "2020 Dec 15", "TIC 262400835, 2020 Dec 15"),
    ("TIC 262400835", "2020 Dec 16", "TIC 262400835, 2020 Dec 16"),
    ("LP 12-502",     "2023 Dec 08", "LP 12-502, 2023 Dec 08"),
    ("LP 12-502",     "2023 Dec 15", "LP 12-502, 2023 Dec 15"),
    ("TIC 300651846", "2025 Feb 10", "TIC 300651846, 2025 Feb 11"),
    ("TIC 435899024", "2021 May 25 - 2021 June 01",
     "TIC 435899024, 2021 May 25–Jun 01"),
    ("PTFO 8-8695",   "2014 Jan 19", "PTFO 8-8695, 2014 Jan 19"),
    ("PTFO 8-8695",   "2014 Feb 23", "PTFO 8-8695, 2014 Feb 23"),
    ("PTFO 8-8695",   "2014 Dec 27", "PTFO 8-8695, 2014 Dec 27"),
    ("PTFO 8-8695",   "2015 Jan 10", "PTFO 8-8695, 2015 Jan 10"),
    ("PTFO 8-8695",   "2015 Jan 23", "PTFO 8-8695, 2015 Jan 23"),
    ("PTFO 8-8695",   "2015 Feb 10", "PTFO 8-8695, 2015 Feb 10"),
    ("PTFO 8-8695",   "2015 Feb 14", "PTFO 8-8695, 2015 Feb 14"),
    ("PTFO 8-8695",   "2015 Feb 23", "PTFO 8-8695, 2015 Feb 23"),
    ("PTFO 8-8695",   "2016 Feb 9",  "PTFO 8-8695, 2016 Feb 9"),
    ("PTFO 8-8695",   "2016 Nov 25", "PTFO 8-8695, 2016 Nov 25"),
    ("PTFO 8-8695",   "2016 Dec 9",  "PTFO 8-8695, 2016 Dec 9"),
]

# Bandpass -> plot color (see module docstring). Keyed by the CSV band string.
BAND_COLOR = {
    "B": "b", "V": "g",
    "g": "g", "g'": "g",
    "r": "r", "r'": "r", "R": "r",
    "i": "darkorange",
    "I": "darkred", "I_C": "darkred", "I+z": "darkred",
    "z": "purple", "z_s": "purple",
    "Tierras": "darkorange",
    "TESS": "k",
    "J": "olive", "H": "sienna",
    "K_s": "darkorange", "NB2.09": "maroon",
}

# Effective wavelength (micron) used only to order the stacked bands (blue->red).
BAND_WL = {
    "B": 0.44, "g": 0.48, "g'": 0.48, "V": 0.55, "r": 0.62, "r'": 0.62,
    "R": 0.65, "i": 0.76, "TESS": 0.79, "I": 0.80, "I_C": 0.80, "I+z": 0.83,
    "Tierras": 0.86, "z": 0.91, "z_s": 0.91, "J": 1.24, "H": 1.63,
    "NB2.09": 2.09, "K_s": 2.15,
}


def band_label(b):
    """TeX-ish label for a band string."""
    special = {
        "TESS": "TESS", "Tierras": "Tierras", "NB2.09": r"2.09 $\mu$m",
        "z_s": r"$z_{\rm s}$", "I_C": r"$I_{\rm C}$", "K_s": r"$K_{\rm s}$",
        "I+z": r"$I\!+\!z$", "g'": r"$g'$", "r'": r"$r'$",
    }
    if b in special:
        return special[b]
    return f"${b}$"


def epoch_bands():
    """Depths CSV -> (groups, refs).

    ``groups`` maps (star, date_ut) -> ordered unique band list; ``refs`` maps
    (star, date_ut) -> the CSV ref number (the data-source citation key).
    """
    rows = list(csv.DictReader(open(export_csv.CSV_PATH)))
    groups, refs = OrderedDict(), {}
    for r in rows:
        if r["star"] in export_csv.DROP_STARS:
            continue
        key = (r["star"], r["date_ut"])
        groups.setdefault(key, [])
        if r["band"] not in groups[key]:
            groups[key].append(r["band"])
        refs.setdefault(key, r["ref"])
    return groups, refs


def load_band(star, date_ut, band):
    """Return the fit + data payload for one (star, date, band) table cell.

    Mirrors fit_dips.export_csv: resolve to a done folded variant if any, and
    prefer the cross-band joint (shared-duration) dips over the single-band
    preferred dips.
    """
    did = export_csv._resolve_id(export_csv._dataset_id(star, date_ut, band))
    rec = registry.load(did)
    name = rec["preferred_model"]
    joint = rec.get("joint_dip_fit")
    if joint and joint.get("dips"):
        dips = joint["dips"]
        base = joint.get("baseline_params", rec["fits"][name]["baseline_params"])
    else:
        dips = rec["fits"][name]["dips"]
        base = rec["fits"][name]["baseline_params"]
    t, flux, _ = loaders.load_lightcurve(rec["data"])
    return {
        "band": band, "t": t, "flux": flux, "name": name, "dips": dips,
        "base": base, "t_ref": rec["t_ref"], "period": rec.get("period_day"),
        "dip_windows": [tuple(w) for w in rec["labels"].get("dips", [])],
        "flare_windows": [tuple(w) for w in rec["labels"].get("flares", [])],
    }


def _in_windows(t, windows):
    m = np.zeros(np.size(t), dtype=bool)
    for x0, x1 in windows:
        m |= (t >= x0) & (t <= x1)
    return m


def model_flux(bd, t, nsigma=0.0, with_dip=True):
    """Model flux for band-dict ``bd`` at times ``t``.

    ``nsigma`` shifts every dip depth by nsigma * its 1-sigma error (clipped at
    >=0), reproducing fit_dips.diagnostics._model_at_depth. ``with_dip=False``
    returns the bare baseline (the no-dip curve).
    """
    _, profile = parse_model_name(bd["name"])
    m = baseline_value(bd["name"], bd["base"], t, bd["t_ref"], bd["period"])
    if not with_dip:
        return m
    flat = []
    for d in bd["dips"]:
        ae = d["A_err"] if np.isfinite(d["A_err"]) else 0.0
        A = max(d["A"] + nsigma * ae, 0.0)
        flat += [A, d["t0"], d["w"], d["tau"]] if profile == "trap" \
            else [A, d["t0"], d["w"]]
    if flat:
        m = m - dip_model(t, flat, profile)
    return m


def timebin(x, y, bw):
    """Mean-bin (x, y) into bins of width ``bw`` in x; empty-safe."""
    if len(x) == 0:
        return np.array([]), np.array([])
    df = pd.DataFrame({"b": np.floor(x / bw), "x": x, "y": y})
    g = df.groupby("b")
    return g["x"].mean().to_numpy(), g["y"].mean().to_numpy()


BINW_HR = 10.0 / 60.0   # 10-min bins for high-cadence bands

# Per-panel x-tick overrides (keyed by panel title).
XTICKS = {
    "LP 12-502, 2023 Dec 08": [-4, -2, 0, 2],
}

# Data-source credit per CSV ``ref`` number (shown as the 2nd title line).
REF_CITE = {
    "1": "This work",
    "2": "Koen (2023)",
    "3": "Yu et al. (2015)",
    "4": "Onitsuka et al. (2017)",
    "5": "Tanimoto et al. (2020)",
}


def draw_panel(ax, title, bands_payload, source=""):
    """Draw one epoch panel: stacked, color-coded per-band overviews."""
    # order bands blue -> red (bluest, deepest on top)
    payloads = sorted(bands_payload, key=lambda b: BAND_WL.get(b["band"], 1.0))

    # common time reference for the panel: mean of all fitted dip mid-times, so
    # simultaneous bands (and the phased TIC 300651846 dips) align at 0 hr.
    all_t0 = [d["t0"] for bd in payloads for d in bd["dips"]]
    t0_ref = float(np.mean(all_t0)) if all_t0 else \
        float(np.median(payloads[0]["t"]))

    # --- first pass: compute each band's robust vertical extent (percent) ---
    spans, xall = [], []
    for bd in payloads:
        t, flux = bd["t"], bd["flux"]
        y = (flux - 1.0) * 100.0
        nonflare = ~_in_windows(t, bd["flare_windows"])
        # grid the model ONLY over the fitted (non-flare) span: the baseline is
        # constrained there and extrapolates wildly across the flare gaps.
        tnf = t[nonflare] if np.any(nonflare) else t
        tg = np.linspace(tnf.min(), tnf.max(), 2000)
        ylo = (model_flux(bd, tg, nsigma=+1) - 1.0) * 100.0
        yhi = (model_flux(bd, tg, nsigma=-1) - 1.0) * 100.0
        ynom = (model_flux(bd, tg) - 1.0) * 100.0
        if np.any(nonflare):
            dlo, dhi = np.nanpercentile(y[nonflare], [1, 99])
        else:
            dlo, dhi = np.nanpercentile(y, [1, 99])
        lo = min(dlo, np.min(ylo), np.min(ynom))
        hi = max(dhi, np.max(yhi), np.max(ynom))
        bd["_prep"] = dict(y=y, nonflare=nonflare, tg=tg, ylo=ylo, yhi=yhi,
                           ynom=ynom, ybase=(model_flux(bd, tg, with_dip=False)
                                             - 1.0) * 100.0, lo=lo, hi=hi)
        spans.append(hi - lo)
        xall.append(((t[nonflare] - t0_ref) * 24.0))
    spacing = max(spans) * 1.12 + 0.4
    if len(payloads) == 2 and title != "PTFO 8-8695, 2014 Jan 19":
        # two-band PTFO epochs: tighten the stack for more per-band dynamic
        # range (2014 Jan 19 excepted -- its wide dips need the old spacing)
        spacing *= 0.75

    # --- second pass: stack and draw ---
    for i, bd in enumerate(payloads):
        p = bd["_prep"]
        off = -spacing * i
        color = BAND_COLOR.get(bd["band"], "0.3")
        t = bd["t"]
        x = (t - t0_ref) * 24.0
        xg = (p["tg"] - t0_ref) * 24.0

        # +/-1 sigma depth band, nominal model, no-dip baseline
        ax.fill_between(xg, p["ylo"] + off, p["yhi"] + off, color=color,
                        alpha=0.28, lw=0, zorder=2)
        ax.plot(xg, p["ynom"] + off, "-", color=color, lw=1.3, zorder=4)
        ax.plot(xg, p["ybase"] + off, ":", color=color, lw=0.9, alpha=0.75,
                zorder=3)

        nf = p["nonflare"]
        y = p["y"]
        dense = bd["band"] in ("TESS", "NB2.09") or t.size > 100
        if dense:
            # gray raw underplot, then black (kept) / gray-x (excluded) binned pts
            ax.scatter(x, y + off, s=1.8, color="0.72", alpha=0.6,
                       linewidths=0, zorder=1, rasterized=True)
            xb, yb = timebin(x[nf], y[nf], BINW_HR)
            ax.scatter(xb, yb + off, s=4.0, color="k", linewidths=0, zorder=6)
            xf, yf = timebin(x[~nf], y[~nf], BINW_HR)
            if len(xf):
                ax.scatter(xf, yf + off, s=11, color="0.45", marker="x",
                           linewidths=0.7, zorder=7)
        else:
            ax.plot(x[nf], y[nf] + off, ".", color="k", ms=3.2, zorder=6)
            ax.plot(x[~nf], y[~nf] + off, "x", color="0.45", ms=3.2, mew=0.8,
                    zorder=7)

        # band label at the right edge (axes-fraction x), near the top of this
        # band's block (data y)
        ax.text(0.985, off + p["hi"] - 0.10 * spacing, band_label(bd["band"]),
                transform=ax.get_yaxis_transform(),
                ha="right", va="top", color=color, fontsize=6.5,
                bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none",
                          alpha=0.75), zorder=1003)

    # axis cosmetics
    xflat = np.concatenate([a for a in xall if len(a)]) if xall else np.array([0])
    xlo, xhi = np.nanmin(xflat), np.nanmax(xflat)
    xpad = 0.06 * (xhi - xlo) + 0.05
    ax.set_xlim(xlo - xpad, xhi + xpad)
    top = payloads[0]["_prep"]["hi"] + 0.10 * spacing
    bot = -spacing * (len(payloads) - 1) + payloads[-1]["_prep"]["lo"] \
        - 0.10 * spacing
    ax.set_ylim(bot, top)
    if title in XTICKS:
        ax.set_xticks(XTICKS[title])
    ax.set_title(f"{title}\n{source}", fontsize=7, pad=3)


def make_figure(panel_specs, groups, refs, outstem):
    set_style("science")
    fig, axs = plt.subplots(3, 3, figsize=(6.0, 9.0))
    axs = axs.ravel()
    for ax, (star, date_ut, title) in zip(axs, panel_specs):
        bands = groups[(star, date_ut)]
        payload = [load_band(star, date_ut, b) for b in bands]
        source = REF_CITE.get(refs.get((star, date_ut)), "")
        draw_panel(ax, title, payload, source)

    # tight_layout in this mpl accounts for sup{x,y}label and reserves a large
    # gap when they sit low, so use plain fig.text() (which it ignores) and
    # place the labels adaptively to hug the measured panel block; the legend is
    # added AFTER layout so it isn't double-counted. bbox_inches='tight' at save
    # trims the outer margins, so only these internal gaps matter.
    fig.tight_layout(pad=0.5, rect=(0.05, 0.09, 1.0, 1.0))
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    H, W = fig.bbox.height, fig.bbox.width
    bbs = [ax.get_tightbbox(r) for ax in axs]
    x0 = min(b.x0 for b in bbs) / W
    y0 = min(b.y0 for b in bbs) / H
    y1 = max(b.y1 for b in bbs) / H
    fig.text(0.53, y0 - 0.012, "Hours from dip center", ha="center", va="top",
             fontsize=13)
    fig.text(x0 - 0.012, 0.5 * (y0 + y1), "Relative flux (%)", rotation=90,
             ha="right", va="center", fontsize=13)

    # shared legend for the marker/line scheme (band colors are labeled in-panel)
    handles = [
        Line2D([], [], color="0.35", lw=1.4, label="preferred model"),
        Patch(facecolor="0.5", alpha=0.35, label=r"depth $\pm1\sigma$"),
        Line2D([], [], color="0.35", lw=1.1, ls=":", label="baseline (no dip)"),
        Line2D([], [], color="k", marker=".", ls="none", ms=6,
               label="in-/out-of-dip"),
        Line2D([], [], color="0.45", marker="x", ls="none", ms=6, mew=1.2,
               label="excluded"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, 0.022), handletextpad=0.4,
               columnspacing=1.3)

    savefig(fig, os.path.join(OUTDIR, f"{outstem}.png"), dpi=350)
    savefig(fig, os.path.join(OUTDIR, f"{outstem}.pdf"), dpi=350)
    plt.close(fig)


def main():
    groups, refs = epoch_bands()
    make_figure(PANELS[:9], groups, refs, "dip_fit_mosaic_1")
    make_figure(PANELS[9:], groups, refs, "dip_fit_mosaic_2")
    print(f"wrote dip_fit_mosaic_1/2 (png+pdf) to {OUTDIR}")


if __name__ == "__main__":
    main()
