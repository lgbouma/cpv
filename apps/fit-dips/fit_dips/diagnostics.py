"""
Diagnostic plots written after a fit, one directory per dataset under
results/fit-dips_output/<dataset_id>/.

For every fitted model:
    <id>__model_<name>.png   two-panel: data+fit (top), residuals (bottom)

Plus dataset-level summaries:
    <id>__overview.png          classification + preferred-model overlay
    <id>__model_comparison.png  dBIC / rescaled & raw reduced-chi2 across models
    <id>__preferred_resid_hist.png  residual histogram for the preferred model

Figures are built with the object-oriented Figure API (not pyplot) so they are
saved straight to disk without ever opening a GUI window.
"""
import glob
import os

import numpy as np
from matplotlib.figure import Figure

from .loaders import REPO_ROOT
from .models import baseline_value, dip_model, parse_model_name

OUTBASE = os.path.join(REPO_ROOT, "results", "fit-dips_output")


def _outdir(dataset_id):
    d = os.path.join(OUTBASE, dataset_id)
    os.makedirs(d, exist_ok=True)
    return d


def _dipflat(mr, nsigma=0.0):
    """Flat dip-parameter list for ``mr``'s profile.

    Each dip depth A is shifted by ``nsigma`` * its 1-sigma error (clipped at
    >= 0); nsigma=0 gives the nominal curve.  Returns (flat_params, profile).
    """
    _, profile = parse_model_name(mr.name)
    out = []
    for d in mr.dips:
        ae = d.A_err if np.isfinite(d.A_err) else 0.0
        A = max(d.A + nsigma * ae, 0.0)
        if profile == "trap":
            out += [A, d.t0, d.w, d.tau]   # [A, t0, W, tau]
        else:
            out += [A, d.t0, d.w]
    return out, profile


def _model_at(name, mr, t, t_ref, period):
    m = baseline_value(name, mr.baseline_params, t, t_ref, period)
    df, profile = _dipflat(mr)
    if df:
        m = m - dip_model(t, df, profile)
    return m


def _baseline_at(name, mr, t, t_ref, period):
    """The baseline model only (no sech dip)."""
    return baseline_value(name, mr.baseline_params, t, t_ref, period)


def _model_at_depth(name, mr, t, t_ref, period, nsigma=0.0):
    """Model with each dip depth shifted by ``nsigma`` * its 1-sigma error.

    Depths are clipped at >= 0 (no emission). nsigma=0 recovers the nominal fit.
    """
    m = baseline_value(name, mr.baseline_params, t, t_ref, period)
    dip, profile = _dipflat(mr, nsigma=nsigma)
    if dip:
        m = m - dip_model(t, dip, profile)
    return m


def _ok(mr):
    return bool(mr.baseline_params) and np.isfinite(mr.BIC)


def _in_dip_mask(t, dip_windows):
    m = np.zeros(np.size(t), dtype=bool)
    for (x0, x1) in dip_windows:
        m |= (t >= x0) & (t <= x1)
    return m


def _save(fig, path):
    fig.savefig(path, dpi=130, bbox_inches="tight")


def _plot_one_model(name, results, t, flux, dip_windows, t_ref, period):
    mr = results[name]
    fig = Figure(figsize=(9, 6))
    ax1, ax2 = fig.subplots(2, 1, sharex=True,
                            gridspec_kw=dict(height_ratios=[3, 1]))

    in_dip = _in_dip_mask(t, dip_windows)
    tgrid = np.linspace(t.min(), t.max(), 1500)
    ok_model = _ok(mr)
    sigma_eff = float(mr.sigma_eff)

    # (1) underplot ALL models beneath the data (thin gray).
    first_other = True
    for nm, other in results.items():
        if nm == name or not _ok(other):
            continue
        ax1.plot(tgrid, _model_at(nm, other, tgrid, t_ref, period), "-",
                 color="0.75", lw=0.8, zorder=1,
                 label="other models" if first_other else None)
        first_other = False

    if ok_model:
        # (2) +/- 1 sigma on the dip depth (semi-transparent).
        ax1.plot(tgrid, _model_at_depth(name, mr, tgrid, t_ref, period, +1), "-",
                 color="tab:green", lw=1.2, alpha=0.4, zorder=2,
                 label=r"depth $\pm1\sigma$")
        ax1.plot(tgrid, _model_at_depth(name, mr, tgrid, t_ref, period, -1), "-",
                 color="tab:green", lw=1.2, alpha=0.4, zorder=2)
        ax1.plot(tgrid, _model_at(name, mr, tgrid, t_ref, period), "-",
                 color="tab:green", lw=1.8, zorder=3, label="this model")
        ax1.plot(tgrid, _baseline_at(name, mr, tgrid, t_ref, period), ":",
                 color="tab:orange", lw=1.5, zorder=3, label="baseline (no dip)")
        for d in mr.dips:
            ax1.axvline(d.t0, color="tab:green", ls=":", lw=1, zorder=2)
        resid = flux - _model_at(name, mr, t, t_ref, period)
        ax2.plot(t, resid, ".", color="k", ms=3, zorder=3)
        ax2.axhspan(-sigma_eff, sigma_eff, color="tab:green", alpha=0.15, lw=0)
    else:
        ax1.text(0.5, 0.5, "model failed to converge", transform=ax1.transAxes,
                 ha="center", va="center", color="tab:red")

    # data on TOP of all model curves.
    ax1.plot(t[~in_dip], flux[~in_dip], ".", color="k", ms=3, zorder=5,
             label="out-of-dip")
    ax1.plot(t[in_dip], flux[in_dip], ".", color="tab:blue", ms=4, zorder=6,
             label="in-dip")

    ax2.axhline(0, color="tab:green", lw=1)
    ax1.legend(loc="best", fontsize=8)
    ax1.set_ylabel("normalized flux")
    ax2.set_ylabel("residual")
    ax2.set_xlabel("time [BJD]")
    ax1.margins(x=0)

    depth_txt = "  ".join(
        f"$\\delta_{i}$={d.delta_pct:.2f}$\\pm${d.delta_pct_err:.2f}% "
        f"FWHM={d.duration*24:.2f}h" for i, d in enumerate(mr.dips))
    ax1.set_title(
        f"{name}   $\\chi^2_\\nu$={mr.chi2_red:.3f} "
        f"($\\chi^2_{{\\nu,0}}$={mr.chi2_red0:.3f})   $\\Delta$BIC={mr.dBIC:.2f}"
        f"   $f$={mr.f:.3f}\n{depth_txt}", fontsize=9)
    return fig


def _plot_overview(dataset_id, t_all, flux_all, dip_windows, flare_windows,
                   best_name, best_mr, t_ref, period, ok_best,
                   model_on_top=False):
    fig = Figure(figsize=(10, 5))
    ax = fig.subplots()

    in_dip = _in_dip_mask(t_all, dip_windows)
    in_flare = _in_dip_mask(t_all, flare_windows)
    in_dip = in_dip & ~in_flare
    oot = ~in_dip & ~in_flare
    for (x0, x1) in dip_windows:
        ax.axvspan(x0, x1, color="tab:blue", alpha=0.08, lw=0)
    for (x0, x1) in flare_windows:
        ax.axvspan(x0, x1, color="tab:red", alpha=0.08, lw=0)

    # By default the data (zorder 5-6) sit above the model curves. For very
    # dense folded light curves (e.g. the full-sector _S88 fold) the model is
    # otherwise buried under the points, so lift the curves above the data.
    z_band, z_model = (7, 9) if model_on_top else (2, 3)

    if ok_best:
        use = ~in_flare
        tgrid = np.linspace(t_all[use].min(), t_all[use].max(), 1500)
        ax.plot(tgrid, _model_at_depth(best_name, best_mr, tgrid, t_ref, period,
                                       +1), "-", color="tab:green", lw=1.2,
                alpha=0.4, zorder=z_band, label=r"depth $\pm1\sigma$")
        ax.plot(tgrid, _model_at_depth(best_name, best_mr, tgrid, t_ref, period,
                                       -1), "-", color="tab:green", lw=1.2,
                alpha=0.4, zorder=z_band)
        ax.plot(tgrid, _model_at(best_name, best_mr, tgrid, t_ref, period), "-",
                color="tab:green", lw=1.8, zorder=z_model,
                label=f"preferred: {best_name}")
        ax.plot(tgrid, _baseline_at(best_name, best_mr, tgrid, t_ref, period),
                ":", color="tab:orange", lw=1.5, zorder=z_model,
                label="baseline (no dip)")

    # data on TOP of the model curves.
    ax.plot(t_all[oot], flux_all[oot], ".", color="k", ms=3, zorder=5,
            label="out-of-dip")
    ax.plot(t_all[in_dip], flux_all[in_dip], ".", color="tab:blue", ms=4,
            zorder=6, label="in-dip")
    ax.plot(t_all[in_flare], flux_all[in_flare], ".", color="tab:red", ms=4,
            zorder=6, label="flare")
    ax.legend(loc="best", fontsize=8)
    ax.set_xlabel("time [BJD]")
    ax.set_ylabel("normalized flux")
    ax.margins(x=0)
    if ok_best and best_mr.dips:
        depth_txt = ", ".join(
            f"$\\delta_{i}$={d.delta_pct:.2f}$\\pm${d.delta_pct_err:.2f}%"
            for i, d in enumerate(best_mr.dips))
    else:
        depth_txt = "no dip"
    ax.set_title(f"{dataset_id}  (preferred = {best_name})   {depth_txt}",
                 fontsize=10)
    return fig


def _plot_comparison(dataset_id, results, preferred):
    names = sorted(results, key=lambda n: results[n].BIC)
    dbic = [results[n].dBIC for n in names]
    chi2 = [results[n].chi2_red for n in names]
    chi2_0 = [results[n].chi2_red0 for n in names]
    colors = ["tab:green" if n == preferred else "0.6" for n in names]
    f = results[preferred].f

    fig = Figure(figsize=(9, 7))
    ax1, ax2, ax3 = fig.subplots(3, 1, sharex=True)
    x = np.arange(len(names))
    ax1.bar(x, dbic, color=colors)
    ax1.set_ylabel(r"$\Delta$BIC")
    ax1.axhline(10, color="k", ls="--", lw=0.8)  # ~decisive threshold
    ax2.bar(x, chi2, color=colors)
    ax2.axhline(1, color="k", ls="--", lw=0.8)
    ax2.set_ylabel(r"$\chi^2_\nu$ (rescaled)")
    ax3.bar(x, chi2_0, color=colors)
    ax3.axhline(1, color="k", ls="--", lw=0.8)
    ax3.set_ylabel(r"$\chi^2_{\nu,0}$ ($\sigma_0$)")
    ax3.set_xticks(x)
    ax3.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
    ax1.set_title(f"{dataset_id}: model comparison (green = preferred; "
                  f"global $f$={f:.3f})", fontsize=10)
    return fig


def _plot_resid_hist(dataset_id, name, mr, t, flux, t_ref, period):
    resid = flux - _model_at(name, mr, t, t_ref, period)
    sigma_eff = float(mr.sigma_eff)
    fig = Figure(figsize=(6, 4))
    ax = fig.subplots()
    ax.hist(resid, bins=max(10, int(np.sqrt(resid.size))), color="0.6",
            density=True)
    xs = np.linspace(resid.min(), resid.max(), 200)
    if sigma_eff > 0:
        ax.plot(xs, np.exp(-0.5 * (xs / sigma_eff) ** 2)
                / (sigma_eff * np.sqrt(2 * np.pi)), "tab:green", lw=1.5,
                label=r"$N(0,\sigma_{\rm eff})$")
    ax.legend(fontsize=8)
    ax.set_xlabel("residual")
    ax.set_ylabel("density")
    ax.set_title(f"{dataset_id}: {name} residuals "
                 f"(rms={np.std(resid):.5f})", fontsize=10)
    return fig


def write_diagnostics(record, t_all, flux_all, dip_windows, flare_windows,
                      t_use, flux_use, t_ref, period, results, preferred):
    """Write all diagnostic figures for a dataset; return the output dir."""
    dataset_id = record["id"]
    outdir = _outdir(dataset_id)
    for old in glob.glob(os.path.join(outdir, "*.png")):
        os.remove(old)

    for name in results:
        fig = _plot_one_model(name, results, t_use, flux_use, dip_windows,
                              t_ref, period)
        _save(fig, os.path.join(outdir, f"{dataset_id}__model_{name}.png"))

    best = results[preferred]
    fig = _plot_overview(dataset_id, t_all, flux_all, dip_windows,
                         flare_windows, preferred, best, t_ref, period,
                         _ok(best), model_on_top=dataset_id.endswith("_S88"))
    _save(fig, os.path.join(outdir, f"{dataset_id}__overview.png"))

    fig = _plot_comparison(dataset_id, results, preferred)
    _save(fig, os.path.join(outdir, f"{dataset_id}__model_comparison.png"))

    if _ok(best):
        fig = _plot_resid_hist(dataset_id, preferred, best, t_use, flux_use,
                               t_ref, period)
        _save(fig, os.path.join(outdir, f"{dataset_id}__preferred_resid_hist.png"))

    return outdir
