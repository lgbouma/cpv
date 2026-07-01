"""Tests for the cross-band shared-duration joint dip fit.

The physical constraint: for a given dip the total duration T14 = 2W and the
mid-time t0 are IDENTICAL across bands, while depth A and ingress/egress
fraction r = tau/W vary per band. These exercise fitting.fit_group_shared_dip
on synthetic trapezoids and the MultiDipLabeler group-fit hook.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pytest  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fit_dips import fitting, models  # noqa: E402
from fit_dips.labeling import MultiDipLabeler  # noqa: E402


def _trap_flux(t, A, t0, W, tau, base=1.0):
    return base - models._trapezoid_dip(t, A, t0, W, tau)


def _band(t, flux, name="poly1_trap", sigma0=1e-3, f=1.0, init=None):
    return {"name": name, "t": t, "flux": flux,
            "dip_windows": [(0.35, 0.65)], "baseline_params": [1.0, 0.0],
            "t_ref": 0.5, "period": None, "sigma0": sigma0, "f": f,
            "init": init}


def test_shared_W_is_identical_across_bands():
    # two bands: SAME W and t0, but different depth and ingress fraction.
    t = np.linspace(0.0, 1.0, 400)
    W, t0 = 0.10, 0.50
    f1 = _trap_flux(t, 0.05, t0, W, 0.02)   # 5%, sharp ingress
    f2 = _trap_flux(t, 0.08, t0, W, 0.05)   # 8%, soft ingress
    res = fitting.fit_group_shared_dip([_band(t, f1), _band(t, f2)])

    assert res is not None and res["converged"]
    sh = res["shared"][0]
    assert sh["W"] == pytest.approx(W, abs=0.01)
    d1 = res["bands"][0]["dips"][0]
    d2 = res["bands"][1]["dips"][0]
    # both bands report the literally-shared W and t0
    assert d1.w == d2.w == pytest.approx(sh["W"])
    assert d1.t0 == d2.t0 == pytest.approx(sh["t0"])
    # depth recovered and band-dependent
    assert d1.delta_pct == pytest.approx(5.0, abs=0.4)
    assert d2.delta_pct == pytest.approx(8.0, abs=0.4)
    # ingress/egress differs per band (free)
    assert d1.tau == pytest.approx(0.02, abs=0.02)
    assert d2.tau == pytest.approx(0.05, abs=0.02)
    assert not np.isclose(d1.tau, d2.tau, atol=1e-3)


def test_presence_mask_skips_band_without_points_in_window():
    # band B only spans [0, 0.30] -> no points under the (0.35, 0.65) window.
    t_a = np.linspace(0.0, 1.0, 400)
    t_b = np.linspace(0.0, 0.30, 60)
    fa = _trap_flux(t_a, 0.05, 0.50, 0.10, 0.02)
    fb = np.ones_like(t_b)
    res = fitting.fit_group_shared_dip([_band(t_a, fa), _band(t_b, fb)])

    assert res is not None
    assert res["bands"][0]["dips"]            # covering band keeps the dip
    assert res["bands"][1]["dips"] == []      # gap band contributes none
    assert res["shared"][0]["W"] == pytest.approx(0.10, abs=0.02)


def test_group_stats_record_initial_and_joint_bic():
    # run the real Stage-1 pipeline so each band carries an init ModelResult,
    # then confirm the joint block reports BOTH stages and a tie dBIC.
    t = np.linspace(0.0, 1.0, 400)
    rng = np.random.default_rng(0)
    f1 = _trap_flux(t, 0.05, 0.50, 0.10, 0.02) + rng.normal(0, 1e-3, t.size)
    f2 = _trap_flux(t, 0.08, 0.50, 0.10, 0.05) + rng.normal(0, 1e-3, t.size)

    bands = []
    for flux in (f1, f2):
        sigma0 = 1e-3
        results, preferred = fitting.fit_all_models(
            t, flux, [(0.35, 0.65)], sigma0, period=None, t_ref=0.5)
        best = results[preferred]
        bands.append(_band(t, flux, name=preferred, sigma0=sigma0, f=best.f,
                           init=best))
    res = fitting.fit_group_shared_dip(bands)

    g = res["group"]
    for key in ("chi2_initial", "chi2_joint", "BIC_initial", "BIC_joint",
                "dBIC_tie", "k_initial", "k_joint"):
        assert key in g and np.isfinite(g[key])
    # the tie removes free parameters, so the joint model has fewer params
    assert g["k_joint"] < g["k_initial"]
    # per-band joint-stage goodness of fit is recorded
    for br in res["bands"]:
        assert np.isfinite(br["chi2_red"]) and np.isfinite(br["BIC"])


def test_group_fit_hook_called_once_and_distributes(monkeypatch):
    t = np.linspace(0.0, 1.0, 51)
    calls = []

    def group_fit(state):
        calls.append(len(state.dips))
        grid = np.linspace(0, 1, 10)
        return [{"tgrid": grid, "model": np.ones_like(grid),
                 "preferred": f"poly1_trap[{i}]", "dip_markers": [0.5],
                 "table_text": f"t{i}"} for i in range(2)]

    panels = [{"t": t, "flux": np.ones_like(t), "title": "a",
               "on_fit": lambda s: None},
              {"t": t, "flux": np.ones_like(t), "title": "b",
               "on_fit": lambda s: None}]
    lab = MultiDipLabeler(panels, "grp", on_group_fit=group_fit)
    lab._on_select(0.4, 0.6)
    lab._classify_preview()
    lab._do_fit()

    assert calls == [1]                       # called exactly once
    assert all(p["model_line"] is not None for p in lab.panels)
    assert lab.state.confirmed
