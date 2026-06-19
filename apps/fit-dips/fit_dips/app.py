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


def build_fit_callback(record):
    """Return on_fit(label_state) that fits, saves, and returns overlay data."""
    t_all, flux_all, _ = loaders.load_lightcurve(record["data"])
    period = record.get("period_day")

    def on_fit(state):
        in_dip, in_flare, oot = state.masks(t_all)
        use = ~in_flare                       # fit over all non-flare points
        t = t_all[use]
        flux = flux_all[use]
        sigma0 = p2p_rms(flux)
        t_ref = float(np.median(t))
        dip_windows = [tuple(w) for w in state.dips]

        results, preferred = fitting.fit_all_models(
            t, flux, dip_windows, sigma0, period=period, t_ref=t_ref)

        # Build overlay curve from the preferred model.
        best = results[preferred]
        tgrid = np.linspace(t.min(), t.max(), 1500)
        _, profile = parse_model_name(preferred)
        dipflat = []
        for d in best.dips:
            dipflat += [d.A, d.t0, d.w, d.tau] if profile == "trap" \
                else [d.A, d.t0, d.w]
        baseline = baseline_value(preferred, best.baseline_params, tgrid, t_ref,
                                  period)
        model = baseline - dip_model(tgrid, dipflat, profile) if dipflat else \
            np.array(baseline)
        dip_markers = [d.t0 for d in best.dips]

        # Diagnostic plots (per-model + summaries).
        outdir = diagnostics.write_diagnostics(
            record, t_all, flux_all, dip_windows, state.flares, t, flux,
            t_ref, period, results, preferred)

        # Persist everything.
        record["sigma0_p2p_rms"] = float(sigma0)
        record["labels"] = {"dips": state.dips, "flares": state.flares,
                            "confirmed": True}
        record["fits"] = {nm: r.to_dict() for nm, r in results.items()}
        record["preferred_model"] = preferred
        record["t_ref"] = t_ref
        record["error_rescale_f"] = float(best.f)
        record["sigma_eff"] = float(best.sigma_eff)
        record["diagnostics_dir"] = outdir
        registry.set_status(record, registry.STATUS_DONE)
        registry.save(record)

        table_text = (_comparison_table(results, preferred, best)
                      + f"\n\ndiagnostics -> {outdir}")
        return {
            "tgrid": tgrid, "model": model, "baseline": baseline,
            "preferred": preferred, "dip_markers": dip_markers,
            "table_text": table_text,
        }

    return on_fit, t_all, flux_all


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


def run_dataset(dataset_id):
    from . import labeling  # imported here so --list works headless
    rec = registry.load(dataset_id)
    if rec.get("data") is None:
        print(f"dataset {dataset_id} has no resolved data file yet "
              f"(note: {rec.get('note','')}).")
        return
    try:
        on_fit, t_all, flux_all = build_fit_callback(rec)
    except Exception as e:  # noqa: BLE001
        registry.set_status(rec, registry.STATUS_FAILED)
        registry.save(rec)
        print(f"failed to load {dataset_id}: {e}")
        raise
    title = (f"{rec['id']}  ({rec['star']} {rec['date_ut']} {rec['band']}, "
             f"{rec['instrument']})")
    labeler = labeling.DipLabeler(t_all, flux_all, title, on_fit)
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


def make_phasefold_dataset(base_id, ncycle=2, fit=True):
    """Seed a +/-ncycle phase-folded variant of ``base_id`` and (optionally) fit it."""
    from . import seed_datasets
    rec = seed_datasets.make_phasefold_variant(base_id, ncycle=ncycle, save=True)
    print(f"created {rec['id']}  (window={rec['data']['window']}, "
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
    ap.add_argument("--list", action="store_true", help="list datasets and exit")
    ap.add_argument("--reset", metavar="ID",
                    help="clear a dataset's fit + diagnostics and mark "
                         "not_started, then exit")
    ap.add_argument("--keep-labels", action="store_true",
                    help="with --reset, keep the dip/flare labels")
    ap.add_argument("--fit-saved", metavar="ID",
                    help="headlessly fit a dataset from its saved labels "
                         "(no GUI), then exit")
    ap.add_argument("--make-phasefold", metavar="BASE_ID",
                    help="create a +/-ncycle phase-folded variant of BASE_ID "
                         "and fit it, then exit")
    ap.add_argument("--ncycle", type=int, default=2,
                    help="cycles per side for --make-phasefold (default 2)")
    args = ap.parse_args(argv)

    if args.list:
        print_list()
        return
    if args.make_phasefold:
        if not registry.exists(args.make_phasefold):
            print(f"no such dataset: {args.make_phasefold}")
            sys.exit(1)
        make_phasefold_dataset(args.make_phasefold, ncycle=args.ncycle)
        return
    if args.fit_saved:
        if not registry.exists(args.fit_saved):
            print(f"no such dataset: {args.fit_saved}")
            sys.exit(1)
        fit_saved(args.fit_saved)
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
    run_dataset(dataset_id)


if __name__ == "__main__":
    main()
