# fit-dips

Interactive tool to measure the **depth, duration, and mid-time** of the
dimming events of complex periodic variables (CPVs), while modeling the smooth
spotted-rotation baseline and excluding flares. Results feed Table 3
(`tab:dipdepths`) of `papers/Bouma_2026_cgcd/`.

Design spec: `docs/superpowers/specs/2026-06-15-fit-dips-design.md`.

## Model

Combined model over the non-flare points:

    F(t) = B(t) - sum_i  D_i(t; profile)

Each labeled dip window -> one dip `D_i`. The dip **profile** is a
**flat-bottomed transit** (symmetric trapezoid, 4 par):

- `trap` : full depth `A_i` for `|t - t0_i| <= W_i - tau_i`, then a linear
  ingress/egress to zero at `|t - t0_i| = W_i`. `W_i` is the half total
  duration, `tau_i` the ingress/egress duration. Nests a box (`tau -> 0`) and a
  triangle (`tau -> W`).

A second profile, `sech` (`A_i * sech((t - t0_i)/w_i)`, the Yu+2015 cusped
profile), still exists in the code but is **no longer fit by default**: across
this sample the sech systematically misfit the (flat-bottomed) dips and inflated
their depths, so the default model comparison is the trapezoid alone. The sech
can still be requested explicitly (e.g. `fit_one_model("poly1", ...)` or by
passing `families=[...]` to `fit_all_models`).

- Reported **depth** = `A_i / B(t0_i)` in percent (the full-depth decrement at
  mid-dip).
- Reported **duration** = FWHM (width at half depth) = `2*W_i - tau_i`
  (trapezoid; `2.6339*w_i` for a sech).
- The (trapezoid) dip is combined with every candidate **baseline** and compared
  by BIC. The default baseline pool is `poly1`, `poly2`, `fourier1`, `fourier2`
  (Fourier at the star's known rotation period; skipped if no period is set);
  `poly3/poly4/fourier3/fourier4` remain available by name but are dropped from
  the default pool because they tend to overfit short single-night windows. So
  the default comparison is the flat-bottomed transit against each local-trend
  baseline. Model names encode the profile: the trapezoid appends `_trap`
  (`poly2_trap`); a bare baseline name (`poly2`) denotes the sech.
- **Decoupled (two-stage) fit:** the baseline is fit to the *out-of-dip* points
  only (exact linear LS), then the dip is fit to the residual. This stops a
  flexible baseline from bowing up under the dip (which would inflate the
  depth); the continuum smoothly bridges the dip.
- Uncertainties use a **single global rescale factor** `f` (no per-model
  jitter): fit all models at `sigma0 = p2p_rms`, then set `sigma_eff = f*sigma0`
  with `f` chosen so the BIC-preferred model has reduced chi^2 = 1. Because the
  preferred model itself depends on `f`, this is solved by a short fixed-point
  iteration (preferred -> f -> preferred). It re-weights goodness-of-fit vs.
  complexity in the BIC and inflates/deflates every reported depth uncertainty
  by `f`. Parameter estimates are unchanged by the rescale; no MCMC.

## Quick start (from this directory)

```bash
# one-time: seed the dataset registry from Table 3.
# Non-destructive: only creates datasets that don't exist yet; existing records
# (and any labels/fits in them) are left untouched.
#   --update  refresh metadata (data path, period, note) but PRESERVE labels/fits
#   --force   reset records to clean -- but records with saved work are skipped
python -m fit_dips.seed_datasets

# see all datasets and their status
python -m fit_dips.app --list

# label + fit one dataset (opens an interactive window)
python -m fit_dips.app --dataset LP12-502_20231215_Tierras

# interactive picker
python -m fit_dips.app

# clear a fit that was marked "done" but is unsatisfactory (re-label from scratch)
python -m fit_dips.app --reset LP12-502_20231215_Tierras
#   add --keep-labels to clear the fit but keep the dip/flare windows
```

## Controls (in the plot window)

| key / action      | effect                                              |
|-------------------|-----------------------------------------------------|
| **drag**          | mark a dip window (blue); each drag = one dip       |
| `x`               | switch to flare mode (drags now mark flares, red)   |
| `d`               | back to dip mode                                     |
| `u`               | undo the last window                                 |
| `r`               | reset all labels (and clear the fit overlay)         |
| `enter`           | classify points (auto out-of-dip) and ask to confirm|
| `y`               | confirm -> fit all baselines, overlay the best one   |
| `n`               | cancel the confirmation, keep editing                |
| `q`               | quit the window                                      |

matplotlib's default key bindings are disconnected in the labeler so they don't
shadow these (e.g. `f`=fullscreen, `r`/`h`=home, `s`=save). Pan/zoom/save/home
still work via the toolbar buttons.

After `y`, the preferred model (solid green) and its baseline-only continuum
(orange dotted) are overlaid, the model-comparison table prints to the terminal,
and the full results (all 8 models + the BIC-preferred one) are saved to
`datasets/<id>.json`; status flips to `done`.

## Diagnostic plots

Each fit also writes figures to
`results/fit-dips_output/<dataset_id>/` (stale PNGs cleared on each run):

- `<id>__model_<name>.png` — two-panel **data+fit (top) / residual (bottom)**,
  one per baseline model, titled with chi2_red / chi2_red0 / dBIC / f / depth.
- `<id>__overview.png` — full light curve with in-dip/flare/out-of-dip
  classification and the preferred-model overlay.
- `<id>__model_comparison.png` — dBIC, rescaled reduced chi2, and raw reduced
  chi2 (at sigma0) across all models (preferred highlighted; global f in title).
- `<id>__preferred_resid_hist.png` — residual histogram vs. a Gaussian for the
  preferred model.

## Generating the manuscript table

```bash
python -m fit_dips.make_table   # writes papers/.../tables/dipdepths.tex
```

One row per fitted dip, from the BIC-preferred model of each `done` dataset.
`ms.tex` should `\input{tables/dipdepths.tex}` (not yet wired — do this once
real fits exist).

## Status of datasets (33 of 33 ready)

All datasets have data on disk:
- TIC 262400835: MuSCAT2 Dec 13 & 15 (g, r, i, z), MuSCAT1 Dec 16 (g, r, z),
  TESS S32 SAP_FLUX windows for Dec 13, 15 & 16.
- LP 12-502: KeplerCam B/g (Dec 08, 15), Tierras Dec 15, TESS S73 windows
  for Dec 08 & 15. (Dec 16 excluded: data quality too poor that night.)
- TIC 300651846: FourStar NB2.09 and simultaneous TESS S88, Feb 10.
- PTFO 8-8695 (Yu+2015 Fig 6, digitized): TRAPPIST I+z, FourStar H. Time is
  hours relative to mid-dip (converted to days). Period left null, so only the
  polynomial baselines run — appropriate for these short single-dip windows.
  Both nights are flat-bottomed and decisively prefer the `trap` profile over
  `sech` (I+z: poly1_trap, δ=2.18±0.07%, ΔBIC 57 over sech; H: poly1_trap,
  δ=1.36±0.02%, ΔBIC 368).

Literature digitizations (all `generic_csv`, period left null → polynomial
baselines only, as for PTFO 8-8695):
- TIC 20178925 (Günther+2022, 2019 Nov 07): SPECULOOS-South r' & z' and
  simultaneous TESS. Time is BJD_TDB − 2458795 (already days), flux is
  band-offset relative flux (`time_scale=1`).
- TIC 435899024 / J1252 (Koen+2023, SAAO 1m, 2021 May 25–Jun 01): B, V, R, I
  phase-folded curves. The phase axis is converted to days via
  `time_scale = P = 0.363355 d`, and the differential-magnitude flux is
  converted to relative flux via the loader's `flux_in_mag` option
  (10**(−0.4·mag)). All four prefer poly1_trap with a clear blue→red depth
  gradient (B 12.44±0.40%, V 12.59±0.49%, R 11.04±0.31%, I 9.17±0.20%).

## Notes / TODO

- Periods (set in the registry): LP 12-502 = 18.5611 hr, TIC 262400835 =
  7.15735 hr, TIC 300651846 = 8.254 hr. Fourier baselines run for all three.
- Data access mirrors the project drivers: MuSCAT via
  `complexrotators.plotting._get_all_tic262400835_data`, FourStar/TESS via
  `drivers/plot_fourstar_vs_tess.py`. FourStar reduction (time +60 min, flux
  -0.76, end-trim) is replicated in `loaders.load_fourstar_xls`.
- TESS night windows are computed at seed time from the simultaneous
  ground-based data (±0.04 d pad).
- `p2p_rms` is vendored in `fit_dips/p2p.py` because `cdips` is not importable
  in this environment.

## Layout

```
fit_dips/
  models.py      sech dip + baseline families + combined model
  fitting.py     least-squares fit, global error rescale f, BIC model comparison
  loaders.py     per-instrument light-curve readers
  registry.py    per-dataset JSON store + status
  seed_datasets.py  seed registry from Table 3
  labeling.py    LabelState (headless) + DipLabeler (matplotlib GUI)
  app.py         picker -> load -> label -> fit -> save
  make_table.py  preferred results -> dipdepths.tex
  p2p.py         vendored p2p_rms
datasets/        one JSON per dataset
tests/           pytest: models, fitting, labeling
```
