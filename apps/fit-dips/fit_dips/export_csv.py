"""
Export the measured dip depths into
``papers/Bouma_2026_cgcd/tables/dipdepths_thiswork_and_literature.csv``.

That CSV is the input to ``drivers/run_beta_fit.py``.  Its non-depth columns
(``star``, ``date_ut``, ``band``, ``ref``, ``refit``, ``category``,
``svo_filter``, ``dip``) are hand-curated metadata -- the display date/band
strings, the SVO filter ids, and the reference numbers are NOT reconstructable
from the registry -- so this exporter treats the existing CSV as the metadata
source of truth and overwrites ONLY ``depth_pct`` and ``depth_err_pct`` in
place, pulling the current value from each row's fit_dips registry dataset.

Row -> dataset mapping is ``f"{star_no_spaces}_{YYYYMMDD}_{band_suffix}"``;
``band_suffix`` resolves the few instrument-specific display labels (TIC 262400835
2020-12-16 shows g'/r'/z_s for MuSCAT1 bands g/r/z; PTFO 8-8695 2016-02-09 g'/r'/z_s
map to the gp/rp/zs ids).  As in ``make_table``, a done full-sector ``_S88`` (then
``_pm2cycle``) fold is preferred over the single-cycle base, and the cross-band
joint (shared-duration) dips are used when present.

Stars in ``DROP_STARS`` are removed entirely.  TIC 20178925 is dropped: its
SSO and TESS data show different dip widths, so the single shared-duration
assumption of the fits is not well represented by the data.

Run:  python -m fit_dips.export_csv [--dry-run]   (from apps/fit-dips/)
"""
import argparse
import csv
import datetime
import os

from . import registry

CSV_PATH = os.path.abspath(os.path.join(
    registry.DATASETS_DIR, "..", "..", "..", "papers", "Bouma_2026_cgcd",
    "tables", "dipdepths_thiswork_and_literature.csv"))

# Stars removed from the table entirely (see module docstring).
DROP_STARS = {"TIC 20178925"}

# Hard floor on the reported depth uncertainty (percent). Fits occasionally
# return a formally tiny (or degenerate, ~0) error; we never report below this.
ERR_FLOOR_PCT = 0.05

# Per-star band-label -> dataset-id-suffix overrides. Bands not listed map to
# themselves. The MuSCAT1 (2020-12-16) g'/r'/z_s labels denote the g/r/z
# datasets; PTFO 8-8695's 2016-02-09 g'/r'/z_s are the gp/rp/zs ids.
_BAND_SUFFIX = {
    "TIC 262400835": {"g'": "g", "r'": "r", "z_s": "z"},
    "PTFO 8-8695": {"g'": "gp", "r'": "rp", "z_s": "zs"},
    "TIC 20178925": {"r'": "rprime", "z'": "zprime"},
}


def _date_ymd(date_ut):
    """CSV display date -> YYYYMMDD (the registry id token). Ranges -> start."""
    s = date_ut.split(" - ")[0].strip().replace("June", "Jun").replace(
        "July", "Jul")
    return datetime.datetime.strptime(s, "%Y %b %d").strftime("%Y%m%d")


def _dataset_id(star, date_ut, band):
    suffix = _BAND_SUFFIX.get(star, {}).get(band, band)
    return f"{star.replace(' ', '')}_{_date_ymd(date_ut)}_{suffix}"


def _resolve_id(dataset_id):
    """Prefer a done folded variant (full-sector _S88, then _pm2cycle)."""
    for suffix in ("_S88", "_pm2cycle"):
        variant = f"{dataset_id}{suffix}"
        if registry.exists(variant) and \
                registry.load(variant).get("status") == registry.STATUS_DONE:
            return variant
    return dataset_id


def _registry_dips(dataset_id):
    """Return the preferred-model (or joint) dips list, or None if not done."""
    rid = _resolve_id(dataset_id)
    if not registry.exists(rid):
        return None
    rec = registry.load(rid)
    pref = rec.get("preferred_model")
    if rec.get("status") != registry.STATUS_DONE or not pref \
            or pref not in rec.get("fits", {}):
        return None
    return (rec.get("joint_dip_fit", {}).get("dips")
            or rec["fits"][pref].get("dips", []))


def export(dry_run=False, err_floor=ERR_FLOOR_PCT):
    """Rewrite the depth columns of the CSV from the registry.

    The reported uncertainty is floored at ``err_floor`` percent.

    Rows are regenerated PER DATASET: each (star, date, band) group in the CSV
    is matched to its registry dataset and re-emitted with **one row per fitted
    dip** -- so a refit that adds (or removes) a dip changes the row count, not
    just the depths. Per-dataset metadata (svo_filter, ref, category, display
    date/band) is cloned from the group's first existing row; ``dip`` is renumbered
    1..N and ``refit`` is set to True for ``Literature`` rows (every literature
    star is now refit in this work) and False for ``This work``.

    Returns (rows_written, rows_dropped, changes) where ``changes`` is a list of
    (id, dip, old_depth, old_err, new_depth, new_err) tuples; ``new_*`` is None
    for an added dip (no prior row) and ``old_*`` is None for a registry dataset
    that is not done (its existing rows are kept untouched).
    """
    with open(CSV_PATH, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    # Group existing rows by dataset, preserving first-appearance order and the
    # original per-dip rows (for the change report).
    order, groups = [], {}
    for r in rows:
        did = _dataset_id(r["star"], r["date_ut"], r["band"])
        if did not in groups:
            order.append(did)
            groups[did] = []
        groups[did].append(r)

    out, changes, dropped = [], [], 0
    for did in order:
        grp = groups[did]
        template = grp[0]
        if template["star"] in DROP_STARS:
            dropped += len(grp)
            continue
        dips = _registry_dips(did)
        if dips is None:
            # Registry not done: keep the existing rows untouched (new == None
            # flags an unresolved row).
            for r in grp:
                changes.append((did, r["dip"], r["depth_pct"],
                                r["depth_err_pct"], None, None))
                out.append(r)
            continue
        old_by_dip = {r["dip"]: r for r in grp}
        for i, d in enumerate(dips):
            dip = str(i + 1)
            new_depth = f"{d['delta_pct']:.2f}"
            new_err = f"{max(d['delta_pct_err'], err_floor):.2f}"
            row = dict(template)
            row["dip"] = dip
            row["depth_pct"], row["depth_err_pct"] = new_depth, new_err
            row["refit"] = str(row.get("category") == "Literature")
            prev = old_by_dip.get(dip)
            old_depth = prev["depth_pct"] if prev else None
            old_err = prev["depth_err_pct"] if prev else None
            if (old_depth, old_err) != (new_depth, new_err):
                changes.append((did, dip, old_depth, old_err, new_depth, new_err))
            out.append(row)

    if not dry_run:
        with open(CSV_PATH, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(out)
    return out, dropped, changes


def _main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true",
                    help="report changes without writing the CSV")
    ap.add_argument("--err-floor", type=float, default=ERR_FLOOR_PCT,
                    help=f"floor on reported depth uncertainty in percent "
                         f"(default {ERR_FLOOR_PCT})")
    args = ap.parse_args()

    out, dropped, changes = export(dry_run=args.dry_run,
                                   err_floor=args.err_floor)
    unresolved = [c for c in changes if c[4] is None]
    updated = [c for c in changes if c[4] is not None]
    for did, dip, od, oe, nd, ne in updated:
        before = "(added)" if od is None else f"{od}+/-{oe}"
        print(f"  {did} dip{dip}: {before} -> {nd}+/-{ne}")
    if unresolved:
        print("\nUNRESOLVED (kept, depth untouched -- registry not done):")
        for did, dip, od, oe, _, _ in unresolved:
            print(f"  {did} dip{dip}: {od}+/-{oe}")
    verb = "would write" if args.dry_run else "wrote"
    print(f"\n{verb} {len(out)} rows ({len(updated)} depths updated, "
          f"{dropped} dropped for {sorted(DROP_STARS)}) to {CSV_PATH}")


if __name__ == "__main__":
    _main()
