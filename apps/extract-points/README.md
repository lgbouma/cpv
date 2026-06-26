# extract-points

Digitize the PTFO 8-8695 light-curve panels from Tanimoto et al. 2020
(PASJ 72, 23), figures 12 and 13, recovering the per-point relative flux vs.
time for the optical **I_C** (filled circles) and IR (**K_s / J / H**, crosses)
series — ignoring the red model curves and the residual subpanels.

Purpose: re-extract Tanimoto's photometry so the dip depths can be refit with
the same procedure used in `drivers/run_beta_fit.py`
(`papers/Bouma_2026_cgcd`, Table 3 /
`dipdepths_thiswork_and_literature.csv`, `ref=6`), rather than relying on their
published depths/uncertainties.

## Run

```bash
/Users/luke/local/miniconda3/envs/cpv/bin/python run_extract.py
```

No third-party CV libraries; uses only `numpy`, `scipy.ndimage`, `PIL`,
`pandas`, `matplotlib` (all in the `cpv` env).

## Outputs (`output/`)

| file | contents |
|------|----------|
| `digitized_points.csv` | all 20 panels, long form |
| `digitized_points_csv_epochs.csv` | only the 9 epochs overlapping Table 3 |
| `panels/<fig>_<epoch>.csv` | one file per panel |
| `compare_f12.png`, `compare_f13.png` | reproduced panel grids from digitized data |
| `compare_side_by_side_table3.png` | original crop vs digitized, the 9 Table-3 epochs |
| `diagnostics/overlay_f12.png`, `overlay_f13.png` | detected points drawn on the originals (QA) |

### CSV columns
`figure, panel_row, panel_col, epoch, bjd_offset, in_csv, band, marker,
bjd_minus_offset, rel_flux, bjd`

## Method

1. **Panel detection** (`panels.py`): spines found as long dark runs; the tall
   flux subplots are the largest vertical gaps. → 12 panels in f12, 8 in f13.
2. **Calibration** (`calibrate.py`): inward major tick marks detected on the
   BJD axis (residual-subplot bottom spine) and the relative-flux axis;
   distinguished from minor ticks by inward length; ticks coinciding with a
   spine recovered by extrapolation. Linear pixel→data fit, R² ≈ 1.0 on every
   panel. BJD tick values from `config.py` (hand-read leftmost value + 0.1
   spacing); flux ticks 1.0, 0.9, ...
3. **Color separation** (`color_masks.py`): black markers = low R,G,B minus a
   dilated red-model mask (kills the JPEG halo).
4. **Cleaning** (`markers.py`): the light-gray legend frame and the bottom-left
   epoch label are detected and masked; a margin clears the axis ticks.
5. **Band extraction** (`markers.py`): the optical (I_C) curve is always the
   upper band and the IR the lower (IR is plotted with a downward display
   offset), so each column is split at its dominant y-gap; single-cluster
   columns are assigned by continuity, and outliers off the rolling-median
   trend are dropped. One (time, flux) point per ~2 px column.

## Caveats

- The IR series (`K_s/J/H`) keeps Tanimoto's **arbitrary downward display
  offset** (so the reproduced panels match the originals); its baseline is
  ~0.92, not 1.0. Dip depth is a *relative* measurement, so the offset does not
  affect a depth refit, but absolute IR flux is offset.
- Markers overlap along dense curves, so points are a faithful per-column
  sampling of the plotted curve rather than one point per original exposure.
- `config.py` holds the hand-read per-panel metadata (epoch, BJD offset, bands,
  leftmost BJD tick). The nine Table-3 epochs are listed in `config.CSV_EPOCHS`.

## Tests

```bash
/Users/luke/local/miniconda3/envs/cpv/bin/python tests/test_pipeline.py
```
Checks panel counts, calibration R², and that extracted flux/time are sane.
