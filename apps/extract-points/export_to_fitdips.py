#!/usr/bin/env python
"""Export red_aware digitized points to fit-dips-ready CSVs.

Writes one CSV per (CSV-epoch, band) for the 9 Tanimoto+2020 epochs that overlap
Table 3 (18 datasets: I_C + one IR band each) into the shared literature-curve
directory the fit-dips app reads from:

  <repo>/data/photometry/literature_extracted/Tanimoto2020_<YYYYMMDD>_<band>.csv

Each CSV has two columns, ``bjd_minus_offset,relative_flux``: time in days within
the night (BJD - the panel's integer offset) and flux on the plotted relative
scale. fit-dips' ``generic_csv`` loader reads these with ``time_scale=1.0`` and,
since there is no error column, estimates the per-point error from the p2p_rms.

These files back the ``PTFO8-8695_<YYYYMMDD>_<band>`` dataset records seeded by
``fit_dips.seed_datasets`` (ref=6, Tanimoto+2020). Written to the stable data dir
-- NOT this app's output/ -- so ``run_extract.py`` (which wipes output/) can't
strand the fit-dips inputs.

We always export the canonical ``red_aware`` extraction.

Usage:  python export_to_fitdips.py
"""
from __future__ import annotations

import os
from datetime import datetime

import pandas as pd

from extract_points.extract import extract_all

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
LIT = os.path.join(REPO, "data", "photometry", "literature_extracted")
METHOD = "red_aware"


def iso_date(epoch: str) -> str:
    """'2014 Feb 23' -> '2014-02-23' (matches the fit-dips date_ut convention)."""
    return datetime.strptime(epoch, "%Y %b %d").strftime("%Y-%m-%d")


def export(method: str = METHOD, outdir: str = LIT):
    """Write per-(epoch, band) CSVs for the Table-3 epochs; return a manifest."""
    os.makedirs(outdir, exist_ok=True)
    written = []
    for r in extract_all(method):
        if not r.in_csv:
            continue
        ymd = iso_date(r.epoch).replace("-", "")
        for bd in r.bands:
            df = pd.DataFrame({
                "bjd_minus_offset": [round(float(x), 6) for x in bd.bjd],
                "relative_flux": [round(float(f), 6) for f in bd.flux],
            })
            path = os.path.join(outdir, f"Tanimoto2020_{ymd}_{bd.band}.csv")
            df.to_csv(path, index=False)
            written.append((r.epoch, bd.band, len(df), os.path.basename(path)))
    return written


if __name__ == "__main__":
    manifest = export()
    print(f"Wrote {len(manifest)} CSVs ({METHOD}) to {LIT}/")
    for epoch, band, n, fn in manifest:
        print(f"  {epoch:12s} {band:4s} {n:4d} pts -> {fn}")
