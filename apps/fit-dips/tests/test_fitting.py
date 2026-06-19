import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fit_dips import models, fitting  # noqa: E402


def _make_dataset(baseline_name, baseline_params, dips, period=None,
                  noise=0.002, n=400, seed=0):
    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(0, 0.25, n))  # ~6 hr night
    t_ref = float(np.median(t))
    flux = models.full_model(baseline_name, baseline_params,
                             np.array(dips).ravel(), t, t_ref, period)
    flux = flux + rng.normal(0, noise, n)
    return t, flux, t_ref


def _make_trap_dataset(baseline_name, baseline_params, A, t0, W, tau,
                       period=None, noise=0.0015, n=500, seed=0):
    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(0, 0.25, n))
    t_ref = float(np.median(t))
    base = models.baseline_value(baseline_name, baseline_params, t, t_ref, period)
    flux = base - models._trapezoid_dip(t, A, t0, W, tau)
    flux = flux + rng.normal(0, noise, n)
    return t, flux, t_ref


def test_trap_preferred_on_flat_bottomed_dip():
    # A genuinely flat-bottomed (sharp-edged) dip: the trapezoid must win the
    # BIC comparison over sech and recover the true depth without inflation.
    # (sech is dropped from the default pool, so request it explicitly here to
    # exercise the head-to-head comparison.)
    A, t0, W, tau = 0.04, 0.12, 0.02, 0.004
    t, flux, t_ref = _make_trap_dataset("poly1", [1.0, 0.05], A, t0, W, tau,
                                        seed=0)
    win = [(t0 - 0.035, t0 + 0.035)]
    results, preferred = fitting.fit_all_models(
        t, flux, win, sigma0=0.0015, t_ref=t_ref,
        families=["poly1", "poly1_trap", "poly2", "poly2_trap"])
    # preferred is a trapezoid model
    assert models.parse_model_name(preferred)[1] == "trap"
    assert results[preferred].dips[0].profile == "trap"
    # trapezoid recovers the true depth; the best sech inflates it.
    trap_depth = results[preferred].dips[0].delta_pct
    best_sech = min((r for nm, r in results.items()
                     if models.parse_model_name(nm)[1] == "sech"),
                    key=lambda r: r.BIC)
    assert abs(trap_depth - 100 * A) < 0.6
    assert best_sech.dips[0].delta_pct > trap_depth + 0.3   # sech overshoots
    # tau / FWHM are populated for the trapezoid; depth error is finite.
    d = results[preferred].dips[0]
    assert d.tau > 0 and np.isfinite(d.duration) and d.delta_pct_err > 0


def test_sech_preferred_on_sech_dip_when_requested():
    # The sech profile still works and wins on a true sech dip when explicitly
    # included (it is no longer in the default pool).
    A, t0, w = 0.04, 0.12, 0.012
    t, flux, t_ref = _make_dataset("poly1", [1.0, 0.05], [[A, t0, w]],
                                   noise=0.0015, n=500, seed=3)
    win = [(t0 - 0.035, t0 + 0.035)]
    results, preferred = fitting.fit_all_models(
        t, flux, win, sigma0=0.0015, t_ref=t_ref,
        families=["poly1", "poly1_trap", "poly2", "poly2_trap"])
    assert models.parse_model_name(preferred)[1] == "sech"
    assert results[preferred].dips[0].tau == 0.0


def test_default_pool_is_trap_only_low_order():
    A, t0, W, tau = 0.03, 0.12, 0.02, 0.004
    t, flux, t_ref = _make_trap_dataset("poly1", [1.0, 0.1], A, t0, W, tau,
                                        seed=3)
    results, preferred = fitting.fit_all_models(
        t, flux, [(t0 - 0.04, t0 + 0.04)], sigma0=0.0015, t_ref=t_ref)
    # period None => no fourier; default drops poly3/poly4 AND sech -> trap only.
    assert set(results) == {"poly1_trap", "poly2_trap"}
    assert models.parse_model_name(preferred)[1] == "trap"


def test_recovers_single_dip_on_linear_baseline():
    A, t0, w = 0.03, 0.12, 0.01
    t, flux, t_ref = _make_dataset("poly1", [1.0, 0.05], [[A, t0, w]],
                                   noise=0.001)
    res = fitting.fit_one_model(
        "poly1", t, flux, [(t0 - 0.04, t0 + 0.04)], sigma0=0.001, t_ref=t_ref)
    assert res.converged
    dip = res.dips[0]
    assert abs(dip.A - A) < 0.005
    assert abs(dip.t0 - t0) < 0.005
    assert abs(dip.w - w) < 0.005
    # depth in percent ~ A / baseline(~1.0) * 100
    assert abs(dip.delta_pct - 100 * A) < 1.0


def test_bic_prefers_simpler_true_baseline():
    # Data from a linear baseline + a (true) trapezoid dip -> the poly1 baseline
    # should win on BIC over poly2 (default pool is trap-only -> "poly1_trap").
    A, t0, W, tau = 0.03, 0.12, 0.02, 0.004
    t, flux, t_ref = _make_trap_dataset("poly1", [1.0, 0.1], A, t0, W, tau,
                                        seed=3)
    results, preferred = fitting.fit_all_models(
        t, flux, [(t0 - 0.04, t0 + 0.04)], sigma0=0.0015, t_ref=t_ref)
    # No fourier (period None); poly1 is the true baseline and should be preferred.
    assert "fourier1_trap" not in results
    assert preferred == "poly1_trap"
    assert results["poly1_trap"].dBIC == 0.0


def test_rescale_recovers_true_noise_and_unit_chi2():
    A, t0, w = 0.02, 0.12, 0.012
    true_noise = 0.004
    t, flux, t_ref = _make_dataset("poly1", [1.0, 0.0], [[A, t0, w]],
                                   noise=true_noise, seed=7)
    # Under-report sigma0; the rescale f should make up the difference.
    sigma0 = 0.001
    res = fitting.fit_one_model(
        "poly1", t, flux, [(t0 - 0.04, t0 + 0.04)], sigma0=sigma0, t_ref=t_ref)
    # f > 1 because the data are noisier than sigma0
    assert res.f > 1
    # sigma_eff = f*sigma0 recovers the true scatter
    assert abs(res.sigma_eff - true_noise) < 0.0015
    # single-model rescale forces reduced chi^2 to unity
    assert abs(res.chi2_red - 1.0) < 1e-6


def test_preferred_model_has_unit_chi2_after_rescale():
    A, t0, w = 0.02, 0.12, 0.012
    t, flux, t_ref = _make_dataset("poly1", [1.0, 0.1], [[A, t0, w]],
                                   noise=0.004, seed=11)
    # sigma0 deliberately too small -> global f should be > 1
    results, preferred = fitting.fit_all_models(
        t, flux, [(t0 - 0.04, t0 + 0.04)], sigma0=0.001, t_ref=t_ref)
    f = results[preferred].f
    assert f > 1
    # one global f shared by all models
    assert len({round(r.f, 9) for r in results.values()}) == 1
    # preferred model's rescaled reduced chi^2 is unity; all dBIC >= 0
    assert abs(results[preferred].chi2_red - 1.0) < 1e-3
    assert min(r.dBIC for r in results.values()) == 0.0
    assert all(r.dBIC >= -1e-9 for r in results.values())


def test_single_model_rescale_is_sigma0_independent():
    # For one model, sigma_eff = f*sigma0 = residual RMS, so the renormalized
    # depth error does not depend on the input sigma0 (self-calibrating).
    A, t0, w = 0.02, 0.12, 0.012
    t, flux, t_ref = _make_dataset("poly1", [1.0, 0.0], [[A, t0, w]],
                                   noise=0.004, seed=5)
    win = [(t0 - 0.04, t0 + 0.04)]
    a = fitting.fit_one_model("poly1", t, flux, win, sigma0=0.004, t_ref=t_ref)
    b = fitting.fit_one_model("poly1", t, flux, win, sigma0=0.001, t_ref=t_ref)
    assert np.isclose(a.sigma_eff, b.sigma_eff, rtol=1e-6)
    assert np.isclose(a.dips[0].delta_pct_err, b.dips[0].delta_pct_err, rtol=1e-6)
    # and that sigma_eff equals the residual standard error
    assert a.f * a.sigma0 == a.sigma_eff
    assert b.f > a.f  # smaller sigma0 -> larger f, same product


def test_decoupled_baseline_not_inflated_under_dip():
    # Linear continuum + a clear dip. A flexible poly4 must NOT bow up under the
    # dip (decoupled fit): its baseline at t0 should match the true continuum,
    # so the depth is not inflated.
    true_b = [1.0, 0.2]
    A, t0, w = 0.06, 0.13, 0.012
    t, flux, t_ref = _make_dataset("poly1", true_b, [[A, t0, w]], noise=0.001,
                                   n=400, seed=2)
    res = fitting.fit_one_model("poly4", t, flux, [(t0 - 0.03, t0 + 0.03)],
                                sigma0=0.001, t_ref=t_ref)
    base_t0 = models.baseline_value("poly4", res.baseline_params, [t0], t_ref)[0]
    true_base = models.baseline_value("poly1", true_b, [t0], t_ref)[0]
    assert abs(base_t0 - true_base) < 0.01            # baseline not spiking
    assert abs(res.dips[0].delta_pct - 100 * A / true_base) < 1.0  # depth honest


def test_fourier_runs_with_period():
    period = 0.18
    t, flux, t_ref = _make_dataset(
        "fourier1", [1.0, 0.03, 0.02], [[0.02, 0.09, 0.01]],
        period=period, noise=0.001, seed=1)
    results, preferred = fitting.fit_all_models(
        t, flux, [(0.05, 0.13)], sigma0=0.001, period=period, t_ref=t_ref)
    # default pool is trap-only; fourier baselines run when a period is set.
    assert "fourier1_trap" in results
    assert results["fourier1_trap"].converged
