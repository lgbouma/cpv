# label-coincidences

Interactively mark x-locations in the DBSP Hα-EW + photometry mosaics where a
feature **coincides in both panels** (top: Hα EW, bottom: Δflux), and share
those labels with the plotting drivers so they render matching triangles:

- **dip** — local minima coincide → C0 (blue), alpha 0.3
- **flare** — local maxima coincide → C1 (orange), alpha 0.5

A marker is drawn just inside each panel's edge nearest the EW/flux seam: `^`
10% above the bottom axis of the top (EW) panel, and `v` 10% below the top axis
of the bottom (flux) panel. The *shape* encodes the row, the *color* the kind.

## Run the labeler

```bash
/Users/luke/local/miniconda3/envs/cpv/bin/python run_label.py            # both figures
/Users/luke/local/miniconda3/envs/cpv/bin/python run_label.py --figure mosaic
/Users/luke/local/miniconda3/envs/cpv/bin/python run_label.py --figure supp
```

One interactive window opens per figure (the same columns as the drivers, using
the drivers' own loaders so what you see matches the final figure).

| action | effect |
|--------|--------|
| **left-click** in a column | add a **dip** (blue) at the nearest EW datapoint x |
| **right-click** in a column | add a **flare** (orange) at the nearest EW datapoint x |
| **click on/near a marker** | remove it |
| **`u`** | undo the last add/remove |
| **`w`** | write labels now (also auto-saved after every change) |
| **`q`** / close window | finish this figure and advance |

Clicks snap to the nearest Hα-EW datapoint x, and every change is auto-saved to
`output/coincidence_labels.csv`. Re-running loads existing labels so you can
edit incrementally.

## Regenerate the figures with the labels

```bash
cd ../../drivers
/Users/luke/local/miniconda3/envs/cpv/bin/python plot_dbsp_halpha_mosaic.py
/Users/luke/local/miniconda3/envs/cpv/bin/python plot_dbsp_halpha_mosaic_supp.py
```

Both drivers read `output/coincidence_labels.csv` and draw the ticks. With no
labels file they render exactly as before (no markers).

## Labels file (`output/coincidence_labels.csv`)

Long form, one row per marker:

| column | meaning |
|--------|---------|
| `figure` | `mosaic` or `supp` |
| `colkey` | content-based column id, e.g. `LP_12-502_20241104` (see `labels.colkey`) |
| `kind` | `dip` or `flare` |
| `x` | location in the column's native x-units (hours for time columns, phase for the phase column) |

## Layout

| file | role |
|------|------|
| `run_label.py` | entry point |
| `label_coincidences/labels.py` | CSV persistence + `colkey` (imported by the drivers) |
| `label_coincidences/render.py` | `draw_coincidence_markers` (imported by the drivers) — single source of marker style; call once per axis with `row='top'`/`'bottom'` |
| `label_coincidences/data.py` | builds plot-ready series by reusing the drivers' loaders |
| `label_coincidences/interactive.py` | the matplotlib event-driven labeler |

`labels.py` and `render.py` are the only two modules the drivers import (via a
`sys.path` insert), so the column identity and marker appearance stay in one
place.

## Tests

```bash
/Users/luke/local/miniconda3/envs/cpv/bin/python tests/test_labels.py
```
Checks `colkey` stability, the save/load roundtrip, and that a headless
`LabelSession` builds every column and draws markers.
