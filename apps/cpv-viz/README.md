# cpv-viz

A lightweight visualization portal for complex periodic variables. This initial
landing page renders a sortable, filterable table from the January 2026
literature compilation.

The displayed sectors and N_sectors columns are sourced from the concat table
in `results/tables/jan2026_compilation`. Refresh them periodically with
`/Users/luke/Dropbox/proj/cpv/apps/cpv-viz/scripts/update_tess_sectors.py`
before rebuilding the JSON asset.

Example refresh command:

```bash
/Users/luke/local/miniconda3/envs/cpv/bin/python /Users/luke/Dropbox/proj/cpv/apps/cpv-viz/scripts/update_tess_sectors.py --backup
```

## Weekly pipeline

Run the full refresh (sectors update + JSON rebuild) with:

```bash
/Users/luke/local/miniconda3/envs/cpv/bin/python /Users/luke/Dropbox/proj/cpv/apps/cpv-viz/scripts/run_refresh_pipeline.py --backup
```

The pipeline writes a status file at
`/Users/luke/Dropbox/proj/cpv/apps/cpv-viz/data/refresh_status.json`
with timestamps and per-step status updates.
By default it also refreshes the literature compilation table that powers the
web view; pass `--skip-lit-update` to skip that step. The updater now masks out
future sectors based on the current TESS midtime.

## Data refresh

Use the cpv conda environment to rebuild the JSON asset from the master CSV:

```bash
/Users/luke/local/miniconda3/envs/cpv/bin/python /Users/luke/Dropbox/proj/cpv/apps/cpv-viz/scripts/build_data.py
```

The script writes a date-stamped JSON file plus
`apps/cpv-viz/data/cpv_lit_subset_manifest.json`, which the web app reads to
find the latest asset and display the last-updated date.
It also updates `apps/cpv-viz/data/cpv_lit_subset.json` as a stable fallback.

## PDF manifest

To enable per-star PDF browsing, expose the CPV vetter PDFs via the app and
build a manifest:

```bash
ln -s /Users/luke/local/complexrotators/tars_jan2026_knownCPVs /Users/luke/Dropbox/proj/cpv/apps/cpv-viz/pdfs
/Users/luke/local/miniconda3/envs/cpv/bin/python /Users/luke/Dropbox/proj/cpv/apps/cpv-viz/scripts/build_pdf_manifest.py --base-url ./pdfs
```

If the PDFs live under `LOCALDIR/cpv_finding/tars_jan2026_knownCPVs`, pass that
directory explicitly with `--pdf-dir`.

## Local preview

Serve the repository root and open the app in a browser:

```bash
/Users/luke/local/miniconda3/envs/cpv/bin/python -m http.server 8000
```

Then visit `http://localhost:8000/apps/cpv-viz/`.
