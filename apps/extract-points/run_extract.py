#!/usr/bin/env python
"""Digitize the Tanimoto et al. 2020 PTFO 8-8695 light-curve panels.

Runs the full pipeline with BOTH extraction methods and writes results into
separate subdirectories under output/ so they can be compared:

  output/per_marker/      Option 1: matched-filter per-marker detection
  output/cadence_binned/  Option 2: cadence-matched binning (median y)
  output/red_aware/       Option 3: per-marker after restoring red-occluded interiors

Each subdir contains:
  - digitized_points.csv         all panels, long form
  - digitized_points_csv_epochs.csv  only epochs overlapping Table 3
  - panel_metadata.csv           per-panel axis limits, box bounds, calib R2
  - panels/<fig>_<epoch>.csv     one file per panel
  - compare_f12.png, compare_f13.png   reproduced panels from digitized data
  - compare_side_by_side_table3.png    original vs digitized, Table-3 epochs
  - diagnostics/overlay_f12.png, overlay_f13.png   QA overlays on originals

Plus, at the output/ root:
  - compare_methods_table3.png   original vs both methods, Table-3 epochs

Usage:  python run_extract.py
"""
from __future__ import annotations

import os
import shutil
import pandas as pd

from extract_points.extract import (
    extract_all, to_dataframe, panel_metadata_dataframe, METHODS)
from extract_points.replot import (
    plot_comparison, plot_diagnostic_overlays, plot_side_by_side,
    plot_method_comparison)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")

# method name -> output subdirectory
SUBDIRS = {"per_marker": "per_marker", "cadence": "cadence_binned",
           "red_aware": "red_aware"}


def run_method(method: str, outdir: str):
    """Extract every panel with `method`, write all CSVs and figures to outdir."""
    paneldir = os.path.join(outdir, "panels")
    diag = os.path.join(outdir, "diagnostics")
    os.makedirs(paneldir, exist_ok=True)
    os.makedirs(diag, exist_ok=True)

    results = extract_all(method)
    df = to_dataframe(results)

    df.to_csv(os.path.join(outdir, "digitized_points.csv"), index=False)
    df[df.in_csv].to_csv(
        os.path.join(outdir, "digitized_points_csv_epochs.csv"), index=False)
    panel_metadata_dataframe(results).to_csv(
        os.path.join(outdir, "panel_metadata.csv"), index=False)

    for r in results:
        sub = df[(df.figure == r.figure) & (df.panel_row == r.row) &
                 (df.panel_col == r.col)]
        tag = r.epoch.replace(" ", "_")
        sub.to_csv(os.path.join(paneldir, f"{r.figure}_{tag}.csv"), index=False)

    plot_comparison("f12", results, os.path.join(outdir, "compare_f12.png"))
    plot_comparison("f13", results, os.path.join(outdir, "compare_f13.png"))
    plot_side_by_side(results,
                      os.path.join(outdir, "compare_side_by_side_table3.png"))
    plot_diagnostic_overlays("f12", results, diag)
    plot_diagnostic_overlays("f13", results, diag)

    print(f"  [{method}] {len(results)} panels, {len(df)} points "
          f"({int(df.in_csv.sum())} in CSV epochs) -> {outdir}/")
    return results, df


def main():
    # Clean prior output (inputs live in the project root, not here).
    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)

    results_by_method = {}
    dfs = {}
    for method, sub in SUBDIRS.items():
        print(f"Extracting all panels with method='{method}' ...")
        res, df = run_method(method, os.path.join(OUT, sub))
        results_by_method[method] = res
        dfs[method] = df

    print("Writing cross-method comparison figure ...")
    plot_method_comparison(results_by_method,
                           os.path.join(OUT, "compare_methods_table3.png"))

    # Per-epoch point-count comparison (Table-3 epochs).
    print("\nTable-3 epoch point counts by method (I_C / IR):")
    base = [r for r in results_by_method["per_marker"] if r.in_csv]
    base.sort(key=lambda r: (r.figure, r.row, r.col))
    for r0 in base:
        cells = []
        for m in SUBDIRS:
            r = next(x for x in results_by_method[m]
                     if (x.figure, x.row, x.col) == (r0.figure, r0.row, r0.col))
            ic = r.bands[0].bjd.size
            ir = r.bands[1].bjd.size if len(r.bands) > 1 else 0
            cells.append(f"{m}={ic}/{ir}")
        print(f"  {r0.epoch:12s}  " + "   ".join(cells))
    print(f"\nWrote outputs under {OUT}/")


if __name__ == "__main__":
    main()
