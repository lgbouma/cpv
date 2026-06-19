"""
Least-squares fitting and model comparison for CPV dips.

Error model: a single global rescale factor (no per-model jitter)
----------------------------------------------------------------
The default per-point uncertainty sigma0 = p2p_rms(flux) is a scalar.  Because
it is constant across points, scaling it by any constant factor f does NOT
change the weighted-least-squares solution for the baseline+dip parameters --
so every model is fit exactly once.

We then renormalize the uncertainties by one global fudge factor f, chosen so
that the reduced chi^2 of the BIC-preferred model equals unity:

    sigma_eff = f * sigma0,    f = sqrt(chi2_red of the preferred model at sigma0)

Since the preferred model is itself chosen by BIC (which depends on f), this is
solved by a short fixed-point iteration:

    1. fit each model once at sigma0; record chi2_0 = sum(r^2)/sigma0^2.
    2. f <- 1; repeat: preferred = argmin_BIC(f); f = sqrt(chi2_red0[preferred])
       until the preferred model stops changing (self-consistent).
    3. report chi2_red, lnL, BIC, dBIC at sigma_eff = f*sigma0 for every model.

Renormalizing inflates (or deflates) every reported parameter uncertainty by f,
so the dip-depth error bars reflect the actual scatter about the preferred
model rather than the raw p2p_rms.  No MCMC.
"""
from dataclasses import dataclass, field, asdict
import numpy as np
from scipy.optimize import least_squares

from . import models


@dataclass
class DipResult:
    A: float
    A_err: float
    t0: float
    t0_err: float
    w: float          # sech width, or (trap) W = half the total duration
    w_err: float
    delta_pct: float
    delta_pct_err: float
    duration: float   # FWHM (days)
    duration_err: float
    profile: str = "sech"     # "sech" or "trap"
    tau: float = 0.0          # (trap) ingress/egress duration; 0 for sech
    tau_err: float = 0.0


@dataclass
class ModelResult:
    name: str
    n: int
    k: int
    chi2_red0: float        # reduced chi^2 at sigma0 (pre-renormalization)
    chi2_red: float         # reduced chi^2 at sigma_eff = f*sigma0
    lnL: float
    BIC: float
    dBIC: float = 0.0
    sigma0: float = 0.0
    f: float = 1.0          # global error-rescale factor (same for all models)
    sigma_eff: float = 0.0  # f * sigma0
    converged: bool = True
    baseline_params: list = field(default_factory=list)
    baseline_param_err: list = field(default_factory=list)
    dips: list = field(default_factory=list)  # list of DipResult

    def to_dict(self):
        return asdict(self)


@dataclass
class _RawFit:
    """Result of the (sigma-independent) least-squares fit for one model."""
    name: str
    converged: bool
    n: int
    k: int
    t_ref: float
    baseline_params: list
    dip_params: list      # flat fit params * ndip (sech: [A,t0,w]; trap: [A,t0,W,r])
    perr_unit: np.ndarray  # parameter errors per unit sigma (sqrt diag (J^T J)^-1)
    chi2_0: float          # sum(resid^2) / sigma0^2
    profile: str = "sech"  # dip profile of this model


def _dip_init_from_windows(t, flux, dip_windows, baseline, profile="sech"):
    """Build initial dip-parameter guesses from labeled dip time-windows.

    ``baseline`` is the (per-point) out-of-dip continuum used to estimate the
    initial depth as local_baseline - min(flux) inside each window.

    For the trapezoid the 4th fit parameter is the ingress *fraction* r in
    (0, 1] (tau = r*W), so the optimizer sees only box bounds; W is the half
    total duration.
    """
    baseline = np.asarray(baseline, dtype=float)
    if baseline.ndim == 0:
        baseline = np.full(t.shape, float(baseline))
    p0, lo, hi = [], [], []
    span = float(np.nanmax(t) - np.nanmin(t))
    for (x0, x1) in dip_windows:
        m = (t >= x0) & (t <= x1)
        if not np.any(m):
            continue
        dw = x1 - x0
        t0 = float(t[m][np.argmin(flux[m])])  # deepest point in window
        local_base = float(np.nanmedian(baseline[m]))
        depth = max(local_base - float(np.nanmin(flux[m])), 1e-4)
        if profile == "trap":
            W = max(dw / 2.0, 1e-4)            # half the labeled window
            p0 += [depth, t0, W, 0.25]         # r=0.25 -> modest ingress
            lo += [0.0, x0 - 0.5 * dw, 1e-5, 1e-4]
            hi += [1.0, x1 + 0.5 * dw, span, 1.0]
        else:
            width = max(dw / 4.0, 1e-4)
            p0 += [depth, t0, width]
            lo += [0.0, x0 - 0.5 * dw, 1e-5]
            hi += [1.0, x1 + 0.5 * dw, span]
    return np.array(p0), np.array(lo), np.array(hi)


def _dip_eval(t, fit_params, profile):
    """Evaluate the dip decrement from fit-space parameters.

    For the trapezoid the fit uses the ingress fraction r; convert r -> tau =
    r*W before evaluating the canonical [A, t0, W, tau] trapezoid.
    """
    if profile == "trap":
        p = np.asarray(fit_params, dtype=float).reshape(-1, 4).copy()
        p[:, 3] = p[:, 3] * p[:, 2]            # r -> tau = r * W
        return models.dip_model(t, p.ravel(), "trap")
    return models.dip_model(t, fit_params, profile)


def fit_raw(name, t, flux, dip_windows, sigma0, period=None, t_ref=None):
    """Two-stage fit of one model; returns a sigma-independent _RawFit.

    Stage 1: fit the baseline to OUT-OF-DIP points only (the baseline is linear
    in its parameters, so this is an exact least squares). The continuum is
    therefore determined by the rotational signal alone and cannot trade off
    against the dip.

    Stage 2: fit the sech dip(s) to the residual (baseline - flux) over all
    non-flare points (the model is baseline MINUS the sech).
    """
    t = np.asarray(t, dtype=float)
    flux = np.asarray(flux, dtype=float)
    if t_ref is None:
        t_ref = float(np.median(t))

    _, profile = models.parse_model_name(name)
    nb = models.n_baseline_params(name)
    npp = models.n_dip_params(profile)
    n = t.size

    in_dip = np.zeros(n, dtype=bool)
    for (x0, x1) in dip_windows:
        in_dip |= (t >= x0) & (t <= x1)
    oot = ~in_dip

    # --- Stage 1: baseline from out-of-dip points (exact linear LS) ----------
    X = models.baseline_design(name, t, t_ref, period)
    Xb, yb = (X[oot], flux[oot]) if np.sum(oot) >= nb else (X, flux)
    b, *_ = np.linalg.lstsq(Xb, yb, rcond=None)
    try:
        bperr_unit = np.sqrt(np.clip(np.diag(np.linalg.inv(Xb.T @ Xb)),
                                     0.0, np.inf))
    except np.linalg.LinAlgError:
        bperr_unit = np.full(nb, np.nan)
    baseline_all = X @ b

    # --- Stage 2: sech dip(s) on the residual (baseline - flux) --------------
    # target is positive inside the dip; model = baseline - sech.
    target = baseline_all - flux
    dip0, dip_lo, dip_hi = _dip_init_from_windows(t, flux, dip_windows,
                                                  baseline_all, profile)
    ndip = dip0.size // npp
    converged = True
    if ndip:
        def dres(dp):
            return _dip_eval(t, dp, profile) - target

        sol = least_squares(dres, dip0, bounds=(dip_lo, dip_hi), method="trf",
                            max_nfev=20000)
        dp = sol.x
        converged = bool(sol.success)
        try:
            JTJ = sol.jac.T @ sol.jac
            dperr_unit = np.sqrt(np.clip(np.diag(np.linalg.inv(JTJ)),
                                         0.0, np.inf))
        except np.linalg.LinAlgError:
            dperr_unit = np.full(npp * ndip, np.nan)
        resid = dres(dp)            # = flux - (baseline - dip)
    else:
        dp = np.array([])
        dperr_unit = np.array([])
        resid = flux - baseline_all

    perr_unit = np.concatenate([bperr_unit, dperr_unit])
    chi2_0 = float(np.sum(resid ** 2) / sigma0 ** 2)
    k = nb + npp * ndip + 1  # +1 for the global noise-scale parameter f

    return _RawFit(name=name, converged=converged, n=n, k=k,
                   t_ref=float(t_ref), baseline_params=[float(x) for x in b],
                   dip_params=[float(x) for x in dp], perr_unit=perr_unit,
                   chi2_0=chi2_0, profile=profile)


def _chi2_red0(raw):
    return raw.chi2_0 / max(raw.n - raw.k, 1)


def _bic_at_f(raw, sigma0, f):
    """BIC of a raw fit when uncertainties are rescaled to f*sigma0."""
    n, k = raw.n, raw.k
    var = (f * sigma0) ** 2
    chi2 = raw.chi2_0 / f ** 2
    lnL = -0.5 * (chi2 + n * np.log(2.0 * np.pi * var))
    return k * np.log(n) - 2.0 * lnL


def _select_f_and_preferred(raws, sigma0, maxit=100):
    """Fixed-point iteration for the global rescale f and preferred model.

    f is set so the preferred model's reduced chi^2 is 1; the preferred model is
    argmin BIC at that f.  Returns (f, preferred_name).
    """
    fams = [nm for nm, r in raws.items()
            if r.converged and np.isfinite(r.chi2_0)]
    if not fams:
        return 1.0, min(raws, key=lambda n: raws[n].chi2_0)

    f, pref = 1.0, None
    for _ in range(maxit):
        bic = {nm: _bic_at_f(raws[nm], sigma0, f) for nm in fams}
        new = min(bic, key=bic.get)
        if new == pref:                      # self-consistent fixed point
            return f, pref
        pref = new
        cr0 = _chi2_red0(raws[pref])
        f = float(np.sqrt(cr0)) if cr0 > 0 else 1.0

    # Did not converge (rare cycle): pick the model with the lowest BIC when
    # judged at its own self-consistent f, then use that f.
    def self_bic(nm):
        cr0 = _chi2_red0(raws[nm])
        fm = float(np.sqrt(cr0)) if cr0 > 0 else 1.0
        return _bic_at_f(raws[nm], sigma0, fm)
    pref = min(fams, key=self_bic)
    cr0 = _chi2_red0(raws[pref])
    return (float(np.sqrt(cr0)) if cr0 > 0 else 1.0), pref


def _finalize(raw, sigma0, f, period):
    """Build a ModelResult from a raw fit at rescaled uncertainty f*sigma0."""
    n, k = raw.n, raw.k
    name = raw.name
    sigma_eff = f * sigma0
    var = sigma_eff ** 2
    chi2 = raw.chi2_0 / f ** 2
    chi2_red0 = _chi2_red0(raw)
    chi2_red = chi2 / max(n - k, 1)
    lnL = -0.5 * (chi2 + n * np.log(2.0 * np.pi * var))
    BIC = k * np.log(n) - 2.0 * lnL

    perr = np.asarray(raw.perr_unit, dtype=float) * sigma_eff
    bp = raw.baseline_params
    nb = len(bp)
    bperr = perr[:nb]
    dperr = perr[nb:]

    profile = raw.profile
    npp = models.n_dip_params(profile)
    dips = []
    dp = raw.dip_params
    for i in range(len(dp) // npp):
        row = dp[npp * i: npp * i + npp]
        erow = dperr[npp * i: npp * i + npp]
        A, t0 = float(row[0]), float(row[1])
        A_e, t0_e = float(erow[0]), float(erow[1])
        if profile == "trap":
            # fit space is [A, t0, W, r]; canonical reported is tau = r*W.
            W, r = float(row[2]), float(row[3])
            W_e, r_e = float(erow[2]), float(erow[3])
            tau = r * W
            tau_e = float(np.hypot(r * W_e, W * r_e))   # diag propagation
            w_val, w_e = W, W_e
            dur = float(models.trap_fwhm(W, tau))       # 2W - tau
            dur_e = float(np.hypot(2.0 * W_e, tau_e))
        else:
            w_val, w_e = float(row[2]), float(erow[2])
            tau, tau_e = 0.0, 0.0
            dur = float(models.w_to_fwhm(w_val))
            dur_e = float(models.FWHM_OVER_W * w_e)
        base = float(models.baseline_value(name, bp, [t0], raw.t_ref, period)[0])
        dips.append(DipResult(
            A=A, A_err=A_e, t0=t0, t0_err=t0_e,
            w=w_val, w_err=w_e,
            delta_pct=float(100.0 * A / base),
            delta_pct_err=float(100.0 * A_e / base),
            duration=dur, duration_err=dur_e,
            profile=profile, tau=tau, tau_err=tau_e,
        ))

    return ModelResult(
        name=name, n=n, k=k, chi2_red0=float(chi2_red0),
        chi2_red=float(chi2_red), lnL=float(lnL), BIC=float(BIC),
        sigma0=float(sigma0), f=float(f), sigma_eff=float(sigma_eff),
        converged=raw.converged,
        baseline_params=[float(x) for x in bp],
        baseline_param_err=[float(x) for x in bperr], dips=dips)


def fit_one_model(name, t, flux, dip_windows, sigma0, period=None, t_ref=None):
    """Fit a single model, rescaling so its own reduced chi^2 is unity."""
    raw = fit_raw(name, t, flux, dip_windows, sigma0, period=period,
                  t_ref=t_ref)
    cr0 = _chi2_red0(raw)
    f = float(np.sqrt(cr0)) if cr0 > 0 else 1.0
    return _finalize(raw, sigma0, f, period)


def fit_all_models(t, flux, dip_windows, sigma0, period=None, t_ref=None,
                   families=None):
    """Fit every baseline family, choose the global rescale f, and rank by BIC.

    Fourier families are skipped when ``period`` is None.

    Returns
    -------
    results : dict name -> ModelResult (dBIC filled in; all share one f)
    preferred : name of the minimum-BIC model (its chi2_red == 1 by construction)
    """
    if families is None:
        families = models.default_model_names()
    if period is None:
        families = [f for f in families
                    if not models.parse_model_name(f)[0].startswith("fourier")]
    if t_ref is None:
        t_ref = float(np.median(np.asarray(t, dtype=float)))

    raws = {}
    for name in families:
        try:
            raws[name] = fit_raw(name, t, flux, dip_windows, sigma0,
                                 period=period, t_ref=t_ref)
        except Exception:  # noqa: BLE001 - record failure, keep going
            raws[name] = _RawFit(name=name, converged=False,
                                 n=int(np.size(t)), k=0, t_ref=t_ref,
                                 baseline_params=[], dip_params=[],
                                 perr_unit=np.array([]), chi2_0=np.inf)

    f, preferred = _select_f_and_preferred(raws, sigma0)

    results = {nm: _finalize(raws[nm], sigma0, f, period) for nm in raws}
    # preferred = the minimum-BIC model (consistent with the fixed point).
    preferred = min(results, key=lambda nm: results[nm].BIC)
    best_bic = results[preferred].BIC
    for r in results.values():
        r.dBIC = float(r.BIC - best_bic)
    return results, preferred
