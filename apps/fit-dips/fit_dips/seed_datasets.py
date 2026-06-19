"""
Seed the dataset registry from Table 3 ("This work" rows) of the manuscript.

All datasets whose reduced light curves are on disk get a resolved ``data``
block and become immediately operable.  Only the Yu+2015 PTFO 8-8695 literature
curves remain data-pending (must be digitized).

TESS night windows are computed at seed time from the simultaneous ground-based
data so each TESS segment isolates the relevant night.

Run:  python -m fit_dips.seed_datasets [--overwrite]   (from apps/fit-dips/)
"""
import numpy as np

from . import registry, loaders

# Rotation periods (days).
P_LP12_502 = 18.5611 / 24.0
P_TIC262 = 7.15735 / 24.0      # TIC 262400835
P_TIC300 = 8.254 / 24.0        # TIC 300651846

MUSCAT = "data/photometry/MUSCAT"
TIERRAS = "data/photometry/TIERRAS"
KEPLERCAM = "data/photometry/KeplerCam"
LIT = "data/photometry/literature_extracted"
FOURSTAR = ("data/photometry/FourStar_RAW_DATA/fourstar_20250210/PHOTOMETRY/"
            "TIC300651846_aperture_photometry.xls")
# Home-relative ("~/...") so they resolve on any machine; loaders expand "~".
TESS_262_FITS = ("~/.lightkurve/cache/mastDownload/TESS/"
                 "tess2020324010417-s0032-0000000262400835-0200-s/"
                 "tess2020324010417-s0032-0000000262400835-0200-s_lc.fits")
TESS_300_PKL = ("~/local/complexrotators/cpv_finding/spoc2min_debug/"
                "300651846_S0088_120sec_spoc2min_cpv_periodsearch.pkl")
TESS_LP_FITS = ("data/photometry/tess/"
                "tess2023341045131-s0073-0000000402980664-0268-s_lc.fits")

M2_FITS = {
    "20201213": f"{MUSCAT}/tic262400835_201213_achromatic_k.fits",
    "20201215": f"{MUSCAT}/tic262400835_201215_achromatic_k.fits",
}
M1_CSV = {
    "g": f"{MUSCAT}/TIC262400835_201216_muscat_g_c2345_r8.csv",
    "r": f"{MUSCAT}/TIC262400835_201216_muscat_r_c2345_r10.csv",
    "z": f"{MUSCAT}/TIC262400835_201216_muscat_z_c2345_r14.csv",
}

TESS_WINDOW_PAD = 0.04  # days of padding on each side of the ground-based night
TESS_MIN_PTS = 20       # below this, a TESS window is treated as no-coverage

# Per-dataset point masks: indices (into the time-sorted light curve) to drop,
# e.g. first-frame outliers. Survives `seed_datasets --update`.
POINT_MASKS = {
    "LP12-502_20231208_B": [0],  # first frame sits ~1.4% high
}


def _ground_window(data_spec, pad=TESS_WINDOW_PAD):
    """[tmin-pad, tmax+pad] of a simultaneous ground-based light curve."""
    t, _, _ = loaders.load_lightcurve(data_spec)
    return [float(t.min() - pad), float(t.max() + pad)]


def _npts(data_spec):
    try:
        t, _, _ = loaders.load_lightcurve(data_spec)
        return int(t.size)
    except Exception:
        return 0


# Datasets whose data file is missing/insufficient; filled dynamically below.
_dynamic_pending = []


def build_resolved():
    """Return list of (id, star, date, band, instrument, ref, data, period, note)."""
    out = []

    # --- TIC 262400835, MuSCAT2, 2020 Dec 13 & 15 (g, r, i, z) ---
    for ts, date in (("20201213", "2020-12-13"), ("20201215", "2020-12-15")):
        for band in ("g", "r", "i", "z"):
            out.append((
                f"TIC262400835_{ts}_{band}", "TIC 262400835", date, band,
                "MuSCAT2", 1,
                {"loader": "muscat2_fits", "path": M2_FITS[ts], "band": band},
                P_TIC262, ""))

    # --- TIC 262400835, MuSCAT1, 2020 Dec 16 (g, r, z) ---
    for band in ("g", "r", "z"):
        out.append((
            f"TIC262400835_20201216_{band}", "TIC 262400835", "2020-12-16",
            band, "MuSCAT1", 1,
            {"loader": "muscat_csv", "path": M1_CSV[band]}, P_TIC262, ""))

    # --- TIC 262400835, TESS night windows (S32) ---
    tess_nights = {
        "20201213": {"loader": "muscat2_fits", "path": M2_FITS["20201213"],
                     "band": "g"},
        "20201215": {"loader": "muscat2_fits", "path": M2_FITS["20201215"],
                     "band": "g"},
        "20201216": {"loader": "muscat_csv", "path": M1_CSV["g"]},
    }
    for ts, date in (("20201213", "2020-12-13"), ("20201215", "2020-12-15"),
                     ("20201216", "2020-12-16")):
        win = _ground_window(tess_nights[ts])
        spec = {"loader": "tess_fits", "path": TESS_262_FITS, "window": win}
        did = f"TIC262400835_{ts}_TESS"
        if _npts(spec) >= TESS_MIN_PTS:
            out.append((did, "TIC 262400835", date, "TESS", "TESS", 1, spec,
                        P_TIC262,
                        "TESS S32 SAP_FLUX, windowed to the ground night."))
        else:
            _dynamic_pending.append((
                did, "TIC 262400835", date, "TESS", "TESS", 1,
                "TESS S32 FITS has no coverage for this night."))

    # --- LP 12-502, KeplerCam B/g (2023 Dec 08, 15) ---
    # (Dec 16 excluded: data quality too poor that night.)
    for ts, date, note in (("20231208", "2023-12-08", ""),
                           ("20231215", "2023-12-15", "")):
        for band in ("B", "g"):
            out.append((
                f"LP12-502_{ts}_{band}", "LP 12-502", date, band, "KeplerCam", 1,
                {"loader": "keplercam_dat",
                 "path": f"{KEPLERCAM}/LP12-502_{ts}_KeplerCam_{band}.dat"},
                P_LP12_502, note))

    # --- LP 12-502, TESS S73 (simultaneous, 2023 Dec 08 & 15) ---
    for ts, date in (("20231208", "2023-12-08"), ("20231215", "2023-12-15")):
        win = _ground_window({
            "loader": "keplercam_dat",
            "path": f"{KEPLERCAM}/LP12-502_{ts}_KeplerCam_g.dat"})
        out.append((
            f"LP12-502_{ts}_TESS", "LP 12-502", date, "TESS", "TESS", 1,
            {"loader": "tess_fits", "path": TESS_LP_FITS, "window": win},
            P_LP12_502, "TESS S73 SAP_FLUX, windowed to the ground night."))

    # --- LP 12-502, Tierras (2023 Dec 15) ---
    out.append((
        "LP12-502_20231215_Tierras", "LP 12-502", "2023-12-15", "Tierras",
        "Tierras", 1,
        {"loader": "tierras_csv",
         "path": f"{TIERRAS}/20231215_TIC402980664_circular_fixed_ap_phot_13.csv"},
        P_LP12_502, ""))

    # --- TIC 300651846, 2025 Feb 10 (FourStar NB2.09 + simultaneous TESS) ---
    fourstar_spec = {"loader": "fourstar_xls", "path": FOURSTAR,
                     "trim_after": 2460717.87}
    out.append((
        "TIC300651846_20250210_NB2.09", "TIC 300651846", "2025-02-10", "NB2.09",
        "FourStar", 1, fourstar_spec, P_TIC300, ""))
    tess300_win = _ground_window(fourstar_spec)
    out.append((
        "TIC300651846_20250210_TESS", "TIC 300651846", "2025-02-10", "TESS",
        "TESS", 1,
        {"loader": "tess_pkl", "path": TESS_300_PKL, "time_offset": 2457000.0,
         "window": tess300_win},
        P_TIC300, "TESS S88 (BTJD->BJD), windowed to the FourStar night."))

    # --- PTFO 8-8695, Yu+2015 Fig 6 (digitized; time = hours rel. to mid-dip) ---
    out.append((
        "PTFO8-8695_20140119_I+z", "PTFO 8-8695", "2014-01-19", "I+z",
        "Yu+2015 (TRAPPIST)", 3,
        {"loader": "generic_csv",
         "path": f"{LIT}/Yu2015_fig6_trappist.csv",
         "time_col": "t_minus_t0_hr", "flux_col": "relative_flux",
         "time_scale": 1.0 / 24.0},
        None, "Digitized from Yu+2015 Fig 6; time relative to mid-dip."))
    out.append((
        "PTFO8-8695_20140119_H", "PTFO 8-8695", "2014-01-19", "H",
        "Yu+2015 (FourStar)", 3,
        {"loader": "generic_csv",
         "path": f"{LIT}/Yu2015_fig6_fourstar.csv",
         "time_col": "t_minus_t0_hr", "flux_col": "relative_flux",
         "time_scale": 1.0 / 24.0},
        None, "Digitized from Yu+2015 Fig 6; time relative to mid-dip."))

    return out


def make_phasefold_variant(base_id, ncycle=2, save=True):
    """Create a phase-folded multi-cycle variant of a fitted base dataset.

    Motivation: a strictly-simultaneous single-cycle light curve (e.g. a TESS
    window matched to one ground-based night) can be too noisy to constrain the
    dip, because the small aperture sees few photons in one cycle. This builds a
    companion dataset ``<base_id>_pm<ncycle>cycle`` that:

      * widens the time window by ``ncycle`` rotation cycles on each side
        (``2*ncycle + 1`` cycles total -- e.g. ncycle=2 -> five cycles);
      * phase-folds all those cycles onto the base cycle using the KNOWN period
        (handled by the loader's ``fold_period``/``fold_ref``); and
      * PROPAGATES the base in-dip / out-of-dip / in-flare labels to every cycle
        by phase -- since folding maps each cycle's points onto the base cycle,
        the base labels (defined in that cycle) apply unchanged to the folded
        times.

    The fit then runs in flux-vs-phase with ~(2*ncycle+1)x the in-dip sampling,
    tightening the dip depth / width / mid-phase relative to the single cycle.

    The base dataset must already carry dip labels, a period, and a data block.
    """
    base = registry.load(base_id)
    P = base.get("period_day")
    if P is None:
        raise ValueError(f"{base_id} has no period_day; cannot phase-fold")
    data = base.get("data")
    if not data:
        raise ValueError(f"{base_id} has no resolved data block")
    labels = base.get("labels") or {}
    dips = [[float(a), float(b)] for (a, b) in labels.get("dips", [])]
    flares = [[float(a), float(b)] for (a, b) in labels.get("flares", [])]
    if not dips:
        raise ValueError(f"{base_id} has no dip labels to propagate")

    # Fold reference = center of the labeled region, so the in-dip and in-flare
    # windows sit symmetrically about phase 0 and stay within +/- P/2 (no wrap).
    edges = [x for w in (dips + flares) for x in w]
    lab_lo, lab_hi = min(edges), max(edges)
    if (lab_hi - lab_lo) >= P:
        raise ValueError(
            f"labeled region spans {(lab_hi - lab_lo) / P:.2f} cycles (>= 1); "
            "phase-folding would wrap the labels onto themselves")
    fold_ref = 0.5 * (lab_lo + lab_hi)

    # Extended window: the base data window (or labeled region) grown by
    # ncycle cycles on each side. The loader keeps only the data that exists.
    win = data.get("window")
    nom_lo, nom_hi = (float(win[0]), float(win[1])) if win else (lab_lo, lab_hi)
    ext = [nom_lo - ncycle * P, nom_hi + ncycle * P]

    new_data = dict(data)
    new_data["window"] = ext
    new_data["fold_period"] = float(P)
    new_data["fold_ref"] = float(fold_ref)

    vid = f"{base_id}_pm{ncycle}cycle"
    rec = registry.new_dataset(
        vid, base["star"], base["date_ut"], base["band"], base["instrument"],
        base["ref"], data=new_data, period_day=P,
        note=(f"Phase-folded +/-{ncycle} cycles of {base_id} "
              f"({2 * ncycle + 1}-cycle window); in-dip/oot/flare labels "
              f"propagated by phase at P={P * 24:.4f} h."))
    rec["labels"] = {"dips": dips, "flares": flares, "confirmed": True}
    if save:
        registry.save(rec)
    return rec


# Datasets still awaiting data (path unresolved); seeded so status is tracked.
PENDING = []


META_FIELDS = ("star", "date_ut", "band", "instrument", "ref", "period_day",
               "note", "data")


def _has_work(rec):
    """True if a record holds user labels or fit results worth preserving."""
    lab = rec.get("labels") or {}
    return bool(lab.get("dips") or lab.get("flares") or rec.get("fits")
                or rec.get("status") == registry.STATUS_DONE)


def _all_entries():
    """Unified (id, star, date, band, inst, ref, data, period, note) list."""
    entries = list(build_resolved())  # populates _dynamic_pending
    for (did, star, date, band, inst, ref, note) in PENDING + _dynamic_pending:
        entries.append((did, star, date, band, inst, ref, None, None, note))
    # Apply per-dataset point masks to the data spec.
    out = []
    for (did, star, date, band, inst, ref, data, period, note) in entries:
        if data is not None and did in POINT_MASKS:
            data = dict(data, mask_indices=list(POINT_MASKS[did]))
        out.append((did, star, date, band, inst, ref, data, period, note))
    return out


def seed(update=False, force=False):
    """Seed/refresh the dataset registry.

    Default (neither flag): only create datasets that don't exist yet; existing
    records are left untouched (NON-destructive).

    update=True: also refresh the metadata fields (data path, period, note, ...)
    of existing records, but PRESERVE labels, fits, status, and preferred model.

    force=True: fully reset existing records to a clean state (DESTRUCTIVE --
    wipes labels and fits). Records with saved work are skipped unless force.
    """
    n_new = n_upd = n_skip = 0
    for (did, star, date, band, inst, ref, data, period, note) in _all_entries():
        if not registry.exists(did):
            rec = registry.new_dataset(did, star, date, band, inst, ref,
                                        data=data, period_day=period, note=note)
            registry.save(rec)
            n_new += 1
            continue
        # Record exists.
        if force:
            if _has_work(registry.load(did)):
                # never silently destroy labeled work, even with --force
                n_skip += 1
                continue
            rec = registry.new_dataset(did, star, date, band, inst, ref,
                                        data=data, period_day=period, note=note)
            registry.save(rec)
            n_upd += 1
        elif update:
            rec = registry.load(did)
            new_meta = dict(star=star, date_ut=date, band=band, instrument=inst,
                            ref=ref, period_day=period, note=note, data=data)
            for k in META_FIELDS:
                rec[k] = new_meta[k]
            registry.save(rec)
            n_upd += 1
        else:
            n_skip += 1
    return n_new, n_upd, n_skip


if __name__ == "__main__":
    import sys
    update = "--update" in sys.argv
    force = "--force" in sys.argv
    n_new, n_upd, n_skip = seed(update=update, force=force)
    print(f"created {n_new}, updated {n_upd}, skipped {n_skip} "
          f"(existing/has-work) in {registry.DATASETS_DIR}")
    print(f"registry now holds {len(registry.list_datasets())} datasets")
    if not (update or force):
        print("note: existing records left untouched. Use --update to refresh "
              "metadata (preserves labels), --force to reset unlabeled records.")
