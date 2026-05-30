# cpv-viz

A lightweight visualization portal for complex periodic variables. This initial
landing page renders a sortable, filterable table from the January 2026
literature compilation.

The displayed sectors and N_sectors columns are sourced from the concat table
in `results/tables/jan2026_compilation`. Refresh them periodically with
`scripts/update_tess_sectors.py` before rebuilding the JSON asset.

## Setup

Set these environment variables once (or add them to your shell profile).
All commands below assume you are running from the `apps/cpv-viz` directory.

```bash
# Python interpreter in the cpv environment
export CPV_PYTHON="${CPV_PYTHON:-python}"
# e.g. export CPV_PYTHON=/path/to/envs/cpv/bin/python

# Source directory for vetter PDFs
export CPV_PDF_SRC="$HOME/local/complexrotators/tars_jan2026_knownCPVs"
```

If the `cpv` conda environment is already activated, `CPV_PYTHON=python` is
sufficient.

## Weekly pipeline

Run the full refresh (sectors update + JSON rebuild) with:

```bash
$CPV_PYTHON scripts/run_refresh_pipeline.py --backup
```

The pipeline writes a status file at `data/refresh_status.json` with
timestamps and per-step status updates. By default it also refreshes the
literature compilation table that powers the web view; pass
`--skip-lit-update` to skip that step. The updater masks out future sectors
based on the current TESS midtime.

## Data refresh

Example sector-only refresh:

```bash
$CPV_PYTHON scripts/update_tess_sectors.py --backup
```

Rebuild the JSON asset from the master CSV:

```bash
$CPV_PYTHON scripts/build_data.py
```

The script writes a date-stamped JSON file plus
`data/cpv_lit_subset_manifest.json`, which the web app reads to find the
latest asset and display the last-updated date. It also updates
`data/cpv_lit_subset.json` as a stable fallback.

## PDF manifest

To enable per-star PDF browsing, expose the CPV vetter PDFs via the app and
build a manifest:

```bash
ln -s "$CPV_PDF_SRC" pdfs
$CPV_PYTHON scripts/build_pdf_manifest.py --base-url ./pdfs
```

If the PDFs live elsewhere, pass that directory explicitly with `--pdf-dir`.

## Local preview

Serve the repository root and open the app in a browser:

```bash
$CPV_PYTHON -m http.server 8000
```

Then visit `http://localhost:8000/apps/cpv-viz/`.
