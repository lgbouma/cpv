"""
Compare the newly re-extracted MuSCAT2 light curves (Murgas 2026-07 delivery,
3-column .dat: BJD, normalized flux, flux_err) against the "old" achromatic
FITS extraction used in apps/fit-dips, for TIC 262400835 on the two nights
2020-12-13 and 2020-12-15.

One figure per night; rows = bands (g, r, i, z_s). Old = grey, new = colored.
Both curves are median-normalized so they overlay on the same scale.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

DATADIR = os.path.expanduser("~/Dropbox/proj/cpv/data/photometry/MUSCAT")
NEWDIR = os.path.join(DATADIR, "TIC262400835_MuSCAT2_Murgas/TIC262400835_MuSCAT2")
OUTDIR = os.path.expanduser("~/Dropbox/proj/cpv/results/muscat2_new_vs_old")
os.makedirs(OUTDIR, exist_ok=True)

# night -> (old fits basename, new .dat date tag)
NIGHTS = {
    "2020-12-13": ("tic262400835_201213_achromatic_k.fits", "20201213"),
    "2020-12-15": ("tic262400835_201215_achromatic_k.fits", "20201215"),
}
# band label -> (old FITS HDU name, new .dat band tag, color)
BANDS = [
    ("g",   "flux_g",   "g",  "tab:blue"),
    ("r",   "flux_r",   "r",  "tab:green"),
    ("i",   "flux_i",   "i",  "tab:orange"),
    ("z_s", "flux_z_s", "zs", "tab:red"),
]


def load_old(fitspath, hdu):
    with fits.open(fitspath) as hl:
        d = hl[hdu].data
        t = np.array(d["time_bjd"], dtype=float)
        f = np.array(d["flux"], dtype=float)
    m = np.isfinite(t) & np.isfinite(f)
    t, f = t[m], f[m]
    f = f / np.nanmedian(f)
    return t, f


def load_new(datpath):
    a = np.genfromtxt(datpath)
    t, f, e = a[:, 0], a[:, 1], a[:, 2]
    m = np.isfinite(t) & np.isfinite(f)
    t, f, e = t[m], f[m], e[m]
    norm = np.nanmedian(f)
    return t, f / norm, e / norm


for night, (oldfits, newtag) in NIGHTS.items():
    fitspath = os.path.join(DATADIR, oldfits)
    fig, axs = plt.subplots(len(BANDS), 1, figsize=(10, 11), sharex=True)

    # common time reference (integer BJD) for the x axis
    t0_all = []
    for _, hdu, btag, _c in BANDS:
        to, _ = load_old(fitspath, hdu)
        t0_all.append(to.min())
    t_ref = int(np.floor(min(t0_all)))

    for ax, (blabel, hdu, btag, color) in zip(axs, BANDS):
        to, fo = load_old(fitspath, hdu)
        datpath = os.path.join(NEWDIR, f"TIC262400835_MuSCAT2_{btag}_{newtag}.dat")
        tn, fn, en = load_new(datpath)

        ax.scatter(to - t_ref, fo, s=9, c="0.55", label=f"old (N={len(to)})",
                   zorder=2, edgecolors="none")
        ax.scatter(tn - t_ref, fn, s=6, c=color, alpha=0.7,
                   label=f"new (N={len(tn)})", zorder=3, edgecolors="none")

        # per-point scatter (p2p rms) for a quick precision comparison
        rms_o = np.nanstd(np.diff(fo)) / np.sqrt(2)
        rms_n = np.nanstd(np.diff(fn)) / np.sqrt(2)
        ax.text(0.01, 0.05,
                f"{blabel}   p2p rms: old={rms_o*1e3:.2f} ppt, new={rms_n*1e3:.2f} ppt",
                transform=ax.transAxes, fontsize=9, va="bottom",
                bbox=dict(fc="white", ec="none", alpha=0.7))
        ax.legend(loc="upper right", fontsize=8, framealpha=0.7, markerscale=2)
        ax.set_ylabel("rel. flux")

    axs[0].set_title(f"TIC 262400835  MuSCAT2  {night}   (old FITS vs new .dat)")
    axs[-1].set_xlabel(f"BJD_TDB $-$ {t_ref}")
    fig.tight_layout()
    outpath = os.path.join(OUTDIR, f"muscat2_new_vs_old_{newtag}.png")
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outpath}")
