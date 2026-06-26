#!/usr/bin/env python
"""Digitize the Tanimoto et al. 2020 PTFO 8-8695 light-curve panels.

Runs the full pipeline and writes, under output/:
  - digitized_points.csv         all panels, long form
  - digitized_points_csv_epochs.csv  only epochs overlapping Table 3
  - panels/<fig>_<epoch>.csv     one file per panel
  - compare_f12.png, compare_f13.png   reproduced panels from digitized data
  - diagnostics/overlay_f12.png, overlay_f13.png   QA overlays on originals

Usage:  python run_extract.py
"""
from __future__ import annotations

import os
import pandas as pd

from extract_points.extract import extract_all, to_dataframe
from extract_points.replot import (
    plot_comparison, plot_diagnostic_overlays, plot_side_by_side)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
PANELDIR = os.path.join(OUT, "panels")
DIAG = os.path.join(OUT, "diagnostics")


def main():
    os.makedirs(PANELDIR, exist_ok=True)
    os.makedirs(DIAG, exist_ok=True)

    print("Extracting all panels ...")
    results = extract_all()
    df = to_dataframe(results)

    df.to_csv(os.path.join(OUT, "digitized_points.csv"), index=False)
    df[df.in_csv].to_csv(
        os.path.join(OUT, "digitized_points_csv_epochs.csv"), index=False)

    for r in results:
        sub = df[(df.figure == r.figure) & (df.panel_row == r.row) &
                 (df.panel_col == r.col)]
        tag = r.epoch.replace(" ", "_")
        sub.to_csv(os.path.join(PANELDIR, f"{r.figure}_{tag}.csv"), index=False)

    n_pts = len(df)
    n_panels = len(results)
    print(f"  {n_panels} panels, {n_pts} digitized points "
          f"({df.in_csv.sum()} in CSV-overlapping epochs)")

    print("Writing comparison figures ...")
    plot_comparison("f12", results, os.path.join(OUT, "compare_f12.png"))
    plot_comparison("f13", results, os.path.join(OUT, "compare_f13.png"))
    plot_side_by_side(results, os.path.join(OUT, "compare_side_by_side_table3.png"))

    print("Writing diagnostic overlays ...")
    plot_diagnostic_overlays("f12", DIAG)
    plot_diagnostic_overlays("f13", DIAG)

    # short per-panel summary
    print("\nPer-panel point counts:")
    for r in results:
        counts = ", ".join(f"{bd.band}:{bd.bjd.size}" for bd in r.bands)
        flag = " [Table 3]" if r.in_csv else ""
        print(f"  {r.figure} {r.epoch:12s} ({counts}){flag}")
    print(f"\nWrote outputs under {OUT}/")


if __name__ == "__main__":
    main()
