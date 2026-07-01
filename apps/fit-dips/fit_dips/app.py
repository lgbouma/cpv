"""
fit-dips: interactive labeling + fitting app for CPV dimming events.

Usage (from apps/fit-dips/):
    python -m fit_dips.app --list
    python -m fit_dips.app --dataset TIC262400835_20201216_g
    python -m fit_dips.app                 # interactive picker
    python -m fit_dips.app --reset LP12-502_20231215_Tierras   # clear a bad fit

Workflow: pick a dataset -> plot flux vs time -> drag to label dips -> 'x' to
label flares -> enter to classify (auto out-of-dip) -> 'y' to fit all baseline
models, compare by BIC, overlay the preferred model, and save results to the
dataset JSON.
"""
import argparse
import os
import shutil
import sys
from dataclasses import asdict

import numpy as np

from . import registry, loaders, fitting, diagnostics
from .p2p import p2p_rms
from .models import baseline_value, dip_model, parse_model_name


def _status_tag(rec):
    has_data = rec.get("data") is not None
    return f"[{rec['status']:13s}]{'' if has_data else ' (no data)'}"


def print_list():
    rows = registry.list_datasets()
    print(f"\n{len(rows)} datasets in registry:\n")
    for i, rec in enumerate(rows):
        per = "" if rec.get("period_day") is None else \
            f" P={rec['period_day']*24:.2f}h"
        print(f"  {i:2d}  {_status_tag(rec)}  {rec['id']:32s} "
              f"{rec['star']:14s} {rec['date_ut']} {rec['band']}{per}")
    print()
    return rows


def pick_interactive():
    rows = print_list()
    sel = input("Enter dataset number or id (blank to quit): ").strip()
    if not sel:
        return None
    if sel.isdigit():
        idx = int(sel)
        if 0 <= idx < len(rows):
            return rows[idx]["id"]
        print("out of range")
        return None
    return sel


def _stage1_fit(record, state):
    """Run the per-band Stage-1 fit (baseline comparison + independent dips).

    Returns a dict bundling the inputs and the model-comparison results; does
    NOT write diagnostics or persist anything, so it can be reused both for the
    single-band path and as the first stage of the cross-band joint fit.
    """
    t_all, flux_all, _ = loaders.load_lightcurve(record["data"])
    period = record.get("period_day")
    _, in_flare, _ = state.masks(t_all)
    use = ~in_flare                            # fit over all non-flare points
    t = t_all[use]
    flux = flux_all[use]
    sigma0 = p2p_rms(flux)
    t_ref = float(np.median(t))
    dip_windows = [tuple(w) for w in state.dips]
    results, preferred = fitting.fit_all_models(
        t, flux, dip_windows, sigma0, period=period, t_ref=t_ref)
    return {"t_all": t_all, "flux_all": flux_all, "t": t, "flux": flux,
            "sigma0": sigma0, "t_ref": t_ref, "period": period,
            "dip_windows": dip_windows, "results": results,
            "preferred": preferred, "best": results[preferred]}


def _overlay_payload(name, baseline_params, dips, t, t_ref, period, table_text):
    """Build the GUI overlay dict for a model from a list of DipResult."""
    _, profile = parse_model_name(name)
    tgrid = np.linspace(t.min(), t.max(), 1500)
    dipflat = []
    for d in dips:
        dipflat += [d.A, d.t0, d.w, d.tau] if profile == "trap" \
            else [d.A, d.t0, d.w]
    baseline = baseline_value(name, baseline_params, tgrid, t_ref, period)
    model = baseline - dip_model(tgrid, dipflat, profile) if dipflat else \
        np.array(baseline)
    return {"tgrid": tgrid, "model": model, "baseline": baseline,
            "preferred": name, "dip_markers": [d.t0 for d in dips],
            "table_text": table_text}


def _persist(record, state, s1, outdir, joint=None):
    """Write the Stage-1 results (and optional joint block) to the registry."""
    best = s1["best"]
    record["sigma0_p2p_rms"] = float(s1["sigma0"])
    record["labels"] = {"dips": state.dips, "flares": state.flares,
                        "confirmed": True}
    record["fits"] = {nm: r.to_dict() for nm, r in s1["results"].items()}
    record["preferred_model"] = s1["preferred"]
    record["t_ref"] = s1["t_ref"]
    record["error_rescale_f"] = float(best.f)
    record["sigma_eff"] = float(best.sigma_eff)
    record["diagnostics_dir"] = outdir
    if joint is not None:
        record["joint_dip_fit"] = joint
    else:
        record.pop("joint_dip_fit", None)   # clear any stale joint on a re-fit
    registry.set_status(record, registry.STATUS_DONE)
    registry.save(record)


def build_fit_callback(record):
    """Return on_fit(label_state) that fits, saves, and returns overlay data."""
    t_all, flux_all, _ = loaders.load_lightcurve(record["data"])

    def on_fit(state):
        s1 = _stage1_fit(record, state)
        outdir = diagnostics.write_diagnostics(
            record, s1["t_all"], s1["flux_all"], s1["dip_windows"],
            state.flares, s1["t"], s1["flux"], s1["t_ref"], s1["period"],
            s1["results"], s1["preferred"])
        _persist(record, state, s1, outdir, joint=None)
        best = s1["best"]
        table_text = (_comparison_table(s1["results"], s1["preferred"], best)
                      + f"\n\ndiagnostics -> {outdir}")
        return _overlay_payload(s1["preferred"], best.baseline_params, best.dips,
                                s1["t"], s1["t_ref"], s1["period"], table_text)

    return on_fit, t_all, flux_all


def _joint_summary_text(joint, bres, name):
    """Terminal summary of the joint shared-duration fit for one band."""
    lines = ["", "JOINT (T14=2W and t0 tied across bands):"]
    for sh in joint["shared"]:
        lines.append(f"  dip {sh['dip']}: T14={sh['T14']*24:.3f}"
                     f"+/-{sh['T14_err']*24:.3f}h  "
                     f"t0={sh['t0']:.5f}+/-{sh['t0_err']:.5f}")
    for d in bres["dips"]:
        lines.append(f"  [{name}] depth={d.delta_pct:.2f}+/-"
                     f"{d.delta_pct_err:.2f}%  FWHM={d.duration*24:.3f}h  "
                     f"tau={d.tau*24:.3f}h")
    lines.append(f"  this band (joint): chi2_red={bres['chi2_red']:.3f}  "
                 f"BIC={bres['BIC']:.2f}")
    g = joint["group"]
    if "dBIC_tie" in g:
        tag = "tie FAVORED" if g["dBIC_tie"] < 0 else "tie disfavored"
        lines.append(f"  group initial->joint: dBIC_tie={g['dBIC_tie']:.2f} "
                     f"({tag})")
    return "\n".join(lines)


def build_group_fit_callback(records):
    """Return on_group_fit(state): Stage-1 every band, then one joint dip refit
    with the total duration (T14=2W) and mid-time t0 tied across bands.

    Bands keep their own per-band baseline + BIC selection (Stage-1); the tie
    only re-estimates the dips. Each band's JSON keeps its full Stage-1 results
    AND gains a ``joint_dip_fit`` block (joint dips + joint-stage chi2/BIC and a
    group initial-vs-joint dBIC). Returns one overlay payload per record, in
    order, so the labeler can redraw every panel.
    """
    member_ids = [rec["id"] for rec in records]

    def on_group_fit(state):
        s1s = [_stage1_fit(rec, state) for rec in records]
        # Only trapezoid preferred models with fitted dips can share a duration.
        trap_idx = [i for i, s1 in enumerate(s1s)
                    if parse_model_name(s1["preferred"])[1] == "trap"
                    and s1["best"].dips]
        joint = None
        if len(trap_idx) >= 2 and state.dips:
            jbands = []
            for i in trap_idx:
                s1, best = s1s[i], s1s[i]["best"]
                jbands.append({
                    "name": s1["preferred"], "t": s1["t"], "flux": s1["flux"],
                    "dip_windows": s1["dip_windows"],
                    "baseline_params": best.baseline_params,
                    "t_ref": s1["t_ref"], "period": s1["period"],
                    "sigma0": s1["sigma0"], "f": best.f, "init": best})
            try:
                joint = fitting.fit_group_shared_dip(jbands)
            except Exception as e:  # noqa: BLE001 - fall back to independent
                print(f"joint shared-duration fit failed ({e}); "
                      f"using independent dips")
                joint = None
        else:
            print("shared-duration tie skipped (need >=2 trapezoid bands "
                  "with dips); fitting bands independently")

        payloads = []
        for i, rec in enumerate(records):
            s1, best = s1s[i], s1s[i]["best"]
            outdir = diagnostics.write_diagnostics(
                rec, s1["t_all"], s1["flux_all"], s1["dip_windows"],
                state.flares, s1["t"], s1["flux"], s1["t_ref"], s1["period"],
                s1["results"], s1["preferred"])
            table_text = _comparison_table(s1["results"], s1["preferred"], best)
            rec_joint, dips_overlay = None, best.dips
            if joint is not None and i in trap_idx:
                bres = joint["bands"][trap_idx.index(i)]
                rec_joint = {
                    "group": member_ids,
                    "preferred_model": s1["preferred"],
                    "baseline_params": [float(x) for x in best.baseline_params],
                    "converged": joint["converged"],
                    "shared": joint["shared"],
                    "dips": [asdict(d) for d in bres["dips"]],
                    "n": bres["n"], "k": bres["k"],
                    "chi2_red0": bres["chi2_red0"], "chi2_red": bres["chi2_red"],
                    "lnL": bres["lnL"], "BIC": bres["BIC"],
                    "sigma0": bres["sigma0"], "f": bres["f"],
                    "sigma_eff": bres["sigma_eff"],
                    "group_stats": joint["group"]}
                dips_overlay = bres["dips"]
                diagnostics.write_joint_diagnostic(rec, s1, best, bres, joint)
                table_text += _joint_summary_text(joint, bres, s1["preferred"])
            table_text += f"\n\ndiagnostics -> {outdir}"
            _persist(rec, state, s1, outdir, joint=rec_joint)
            payloads.append(_overlay_payload(
                s1["preferred"], best.baseline_params, dips_overlay, s1["t"],
                s1["t_ref"], s1["period"], table_text))
        return payloads

    return on_group_fit


def _comparison_table(results, preferred, best):
    lines = [f"sigma0(p2p_rms)={best.sigma0:.5f}  rescale f={best.f:.3f}  "
             f"sigma_eff={best.sigma_eff:.5f}",
             "model      k  chi2_red0  chi2_red    dBIC"]
    for nm in sorted(results, key=lambda n: results[n].BIC):
        r = results[nm]
        mark = " *" if nm == preferred else "  "
        lines.append(f"{nm:9s}{mark}{r.k:3d}  {r.chi2_red0:8.3f}  "
                     f"{r.chi2_red:8.3f}  {r.dBIC:8.2f}")
    lines.append("")
    lines.append(f"preferred = {preferred} (chi2_red=1 by renormalization)")
    for i, d in enumerate(best.dips):
        lines.append(f"  dip {i}: depth={d.delta_pct:.2f}+/-{d.delta_pct_err:.2f}%"
                     f"  t0={d.t0:.5f}  dur(FWHM)={d.duration*24:.3f}h")
    return "\n".join(lines)


def find_epoch_group(primary_id, group=True):
    """Ordered dataset ids that share the primary's epoch and time axis.

    Returns ``[primary_id, *siblings]``: same-``star`` same-``date_ut`` datasets
    that have data AND whose loaded time interval intersects the primary's.
    Two guards keep the group to real, co-temporal observations:
      * phase-folded variants (data carrying a ``fold_ref``, e.g. ``_S88`` /
        ``_pm2cycle``) are skipped -- they are derived products seeded from a
        parent's labels and fit headlessly, never co-labeled. This applies to
        both siblings AND the primary itself: launching a folded variant
        returns just that variant (its folded time axis would otherwise
        spuriously overlap the night's raw windows); and
      * the time-interval intersection drops non-overlapping windows.
    What survives is genuinely simultaneous bands (PTFO I_C/K_s, MuSCAT
    g/r/i/z, an overlapping TESS window). With ``group=False`` returns just the
    primary, so the caller falls back to single-panel labeling.
    """
    prec = registry.load(primary_id)
    if not group or prec.get("data") is None:
        return [primary_id]
    if "fold_ref" in prec["data"]:   # a folded variant is always labeled solo
        return [primary_id]
    p_t, _, _ = loaders.load_lightcurve(prec["data"])
    p_lo, p_hi = float(np.min(p_t)), float(np.max(p_t))
    sibs = []
    for rec in registry.list_datasets():
        if rec["id"] == primary_id or rec.get("data") is None:
            continue
        if (rec.get("star") != prec.get("star")
                or rec.get("date_ut") != prec.get("date_ut")):
            continue
        if "fold_ref" in rec["data"]:   # derived phase-folded variant
            continue
        try:
            s_t, _, _ = loaders.load_lightcurve(rec["data"])
        except Exception:  # noqa: BLE001  -- a sibling we can't load isn't fatal
            continue
        s_lo, s_hi = float(np.min(s_t)), float(np.max(s_t))
        if s_hi >= p_lo and s_lo <= p_hi:   # intervals intersect
            sibs.append(rec["id"])
    return [primary_id] + sorted(sibs)


def run_dataset(dataset_id, group=True):
    from . import labeling  # imported here so --list works headless
    ids = find_epoch_group(dataset_id, group=group)
    if len(ids) > 1:
        print(f"epoch group ({len(ids)} bands): {', '.join(ids)}")

    panels, records, group_title = [], [], None
    for did in ids:
        rec = registry.load(did)
        if rec.get("data") is None:
            print(f"dataset {did} has no resolved data file yet "
                  f"(note: {rec.get('note','')}).")
            if did == dataset_id:
                return
            continue
        try:
            on_fit, t_all, flux_all = build_fit_callback(rec)
        except Exception as e:  # noqa: BLE001
            registry.set_status(rec, registry.STATUS_FAILED)
            registry.save(rec)
            print(f"failed to load {did}: {e}")
            if did == dataset_id:   # primary must load; siblings may be skipped
                raise
            continue
        # Only the dense full-sector folded variant gets the phase-binned overlay.
        bin_period = rec.get("period_day") if did.endswith("_S88") else None
        panels.append({
            "t": t_all, "flux": flux_all, "on_fit": on_fit,
            "bin_period": bin_period,
            "title": f"{rec['band']} ({rec['instrument']})  [{did}]",
        })
        records.append(rec)
        if group_title is None:
            group_title = (f"{rec['star']}  {rec['date_ut']}"
                           if len(ids) > 1 else
                           f"{rec['id']}  ({rec['star']} {rec['date_ut']} "
                           f"{rec['band']}, {rec['instrument']})")
    if not panels:
        return
    # >=2 co-temporal bands -> tie the dip duration (T14) and t0 across them.
    on_group_fit = (build_group_fit_callback(records)
                    if len(records) > 1 else None)
    labeler = labeling.MultiDipLabeler(panels, group_title,
                                       on_group_fit=on_group_fit)
    labeler.run()


def fit_saved(dataset_id):
    """Headlessly fit a dataset from its SAVED labels (no GUI).

    Reproduces the safe re-fit pattern: load the record, build the fit callback,
    replay the stored dip/flare labels through it. This fits all baseline
    families, saves results to the registry, and regenerates diagnostics --
    exactly as the GUI's 'y' confirm would. Intended for datasets whose labels
    are already known (e.g. phase-folded variants seeded from a fitted parent).
    """
    from .labeling import LabelState  # local import keeps --list headless
    rec = registry.load(dataset_id)
    if rec.get("data") is None:
        print(f"dataset {dataset_id} has no resolved data file.")
        return None
    labels = rec.get("labels") or {}
    if not labels.get("dips"):
        print(f"dataset {dataset_id} has no saved dip labels to fit.")
        return None
    on_fit, _t, _f = build_fit_callback(rec)
    state = LabelState()
    state.dips = [list(w) for w in labels.get("dips", [])]
    state.flares = [list(w) for w in labels.get("flares", [])]
    payload = on_fit(state)
    print("\n" + payload["table_text"])
    return payload


def refit_group(primary_id, group=True):
    """Headlessly refit a dataset's epoch GROUP from saved labels (no GUI).

    Reuses the same joint pipeline as the interactive 'y' confirm: each band's
    Stage-1 baseline+BIC selection followed by ONE cross-band joint dip refit
    that ties the total duration (T14=2W) and mid-time t0 across bands. The
    PRIMARY's saved dip/flare windows are used as the shared window set for the
    whole group -- mirroring the interactive grouped labeler, which shares one
    window set across bands so the in/out-of-dip times are identical across
    bandpasses on a night. With ``group=False`` (or a lone dataset) this falls
    back to the single-band saved-label refit (``fit_saved``).

    Intended for scripted bulk refits; the primary should be the band whose
    labels best bracket the dip (its window is imposed on the siblings).
    """
    ids = find_epoch_group(primary_id, group=group)
    records = [registry.load(d) for d in ids]
    records = [r for r in records if r.get("data") is not None]
    if not records:
        print(f"{primary_id}: no loadable data in group")
        return
    if len(records) < 2 or not group:
        fit_saved(records[0]["id"])
        return
    prec = records[0]
    labels = prec.get("labels") or {}
    if not labels.get("dips"):
        print(f"{primary_id}: primary has no saved dip labels; cannot refit")
        return
    from .labeling import LabelState
    state = LabelState()
    state.dips = [list(w) for w in labels.get("dips", [])]
    state.flares = [list(w) for w in labels.get("flares", [])]
    print(f"refit group ({len(records)} bands; shared windows from "
          f"{prec['id']}): {', '.join(r['id'] for r in records)}")
    payloads = build_group_fit_callback(records)(state)
    for rec, payload in zip(records, payloads):
        print(f"\n=== {rec['id']} ===")
        print(payload["table_text"])


def make_phasefold_dataset(base_id, ncycle=2, fit=True, suffix=None,
                           full_window=False):
    """Seed a phase-folded variant of ``base_id`` and (optionally) fit it.

    Default: a +/-ncycle window variant. With ``full_window=True``, fold the
    entire available light curve into ``<base_id>_<suffix>`` instead.
    """
    from . import seed_datasets
    rec = seed_datasets.make_phasefold_variant(
        base_id, ncycle=ncycle, save=True, suffix=suffix,
        full_window=full_window)
    print(f"created {rec['id']}  (window={rec['data'].get('window')}, "
          f"fold_ref={rec['data']['fold_ref']:.6f}, "
          f"P={rec['period_day'] * 24:.4f} h)")
    print(f"  propagated labels: {len(rec['labels']['dips'])} dip(s), "
          f"{len(rec['labels']['flares'])} flare(s)")
    if fit:
        fit_saved(rec["id"])
    return rec["id"]


def reset_dataset(dataset_id, keep_labels=False):
    """Clear a dataset's fit, mark it not_started, and remove its diagnostics."""
    registry.reset(dataset_id, keep_labels=keep_labels)
    diag = os.path.join(diagnostics.OUTBASE, dataset_id)
    if os.path.isdir(diag):
        shutil.rmtree(diag)
    print(f"reset {dataset_id} -> not_started "
          f"(labels {'kept' if keep_labels else 'cleared'}; diagnostics removed)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", help="dataset id to load")
    ap.add_argument("--solo", action="store_true",
                    help="with --dataset, label this band alone (do not "
                         "auto-group its same-night siblings into stacked panels)")
    ap.add_argument("--list", action="store_true", help="list datasets and exit")
    ap.add_argument("--reset", metavar="ID",
                    help="clear a dataset's fit + diagnostics and mark "
                         "not_started, then exit")
    ap.add_argument("--keep-labels", action="store_true",
                    help="with --reset, keep the dip/flare labels")
    ap.add_argument("--fit-saved", metavar="ID",
                    help="headlessly fit a dataset from its saved labels "
                         "(no GUI), then exit")
    ap.add_argument("--refit", metavar="ID",
                    help="headlessly refit ID's epoch GROUP from saved labels "
                         "with the cross-band shared-duration joint fit (no "
                         "GUI), then exit; add --solo to refit ID alone")
    ap.add_argument("--make-phasefold", metavar="BASE_ID",
                    help="create a +/-ncycle phase-folded variant of BASE_ID "
                         "and fit it, then exit")
    ap.add_argument("--ncycle", type=int, default=2,
                    help="cycles per side for --make-phasefold (default 2)")
    ap.add_argument("--full-window", action="store_true",
                    help="with --make-phasefold, fold the ENTIRE light curve "
                         "(drop the window) instead of +/-ncycle; "
                         "requires --suffix")
    ap.add_argument("--suffix", metavar="STR",
                    help="with --make-phasefold, name the variant "
                         "<BASE_ID>_<STR> (e.g. S88) instead of pm<ncycle>cycle")
    args = ap.parse_args(argv)

    if args.list:
        print_list()
        return
    if args.make_phasefold:
        if not registry.exists(args.make_phasefold):
            print(f"no such dataset: {args.make_phasefold}")
            sys.exit(1)
        if args.full_window and not args.suffix:
            print("--full-window requires --suffix")
            sys.exit(1)
        make_phasefold_dataset(args.make_phasefold, ncycle=args.ncycle,
                               suffix=args.suffix,
                               full_window=args.full_window)
        return
    if args.fit_saved:
        if not registry.exists(args.fit_saved):
            print(f"no such dataset: {args.fit_saved}")
            sys.exit(1)
        fit_saved(args.fit_saved)
        return
    if args.refit:
        if not registry.exists(args.refit):
            print(f"no such dataset: {args.refit}")
            sys.exit(1)
        refit_group(args.refit, group=not args.solo)
        return
    if args.reset:
        if not registry.exists(args.reset):
            print(f"no such dataset: {args.reset}")
            sys.exit(1)
        reset_dataset(args.reset, keep_labels=args.keep_labels)
        return
    dataset_id = args.dataset or pick_interactive()
    if not dataset_id:
        return
    if not registry.exists(dataset_id):
        print(f"no such dataset: {dataset_id}")
        sys.exit(1)
    run_dataset(dataset_id, group=not args.solo)


if __name__ == "__main__":
    main()
