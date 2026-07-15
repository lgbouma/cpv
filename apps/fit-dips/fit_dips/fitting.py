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

# Trapezoid geometry cap tied to the labeled in-dip window. The dip half-duration
# W is bounded at TRAP_W_CAP_FACTOR * (window half-width), so the total duration
# T14 = 2W stays within this factor of the box the user drew; the mid-time t0 is
# pinned to the window edges. Together these keep a trapezoid from widening or
# drifting out of its label to merge with an adjacent dip. k=1.2 leaves a little
# slack over an exact bracket (k=1). Applies to both the single-band fit and the
# cross-band joint (shared-duration) fit.
TRAP_W_CAP_FACTOR = 1.1


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
            # Cap W at k*(dw/2) and pin t0 to the window, so the trapezoid stays
            # within the labeled box instead of widening/drifting to merge with
            # an adjacent dip (see TRAP_W_CAP_FACTOR).
            W_cap = max(TRAP_W_CAP_FACTOR * 0.5 * dw, 1e-4)
            p0 += [depth, t0, W, 0.25]         # r=0.25 -> modest ingress
            lo += [0.0, x0, 1e-5, 1e-4]
            hi += [1.0, x1, W_cap, 1.0]
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


def fit_group_shared_dip(bands):
    """Joint cross-band dip refit with a SHARED total duration and mid-time.

    Given several co-temporal bands of the SAME night (each already through the
    per-band Stage-1 model comparison so its baseline is fixed), refit all dips
    at once under the physical constraint that, for each dip, the total duration
    ``T14 = 2W`` and the mid-time ``t0`` are IDENTICAL across bands, while the
    depth ``A`` and the ingress/egress fraction ``r = tau/W`` stay free per band.

    This relies on the decoupled (two-stage) fit: each band's baseline was fit to
    out-of-dip points only and is therefore independent of the dip, so tying the
    dip never disturbs any continuum.  Only the dip stage becomes joint.

    Parameters
    ----------
    bands : list of dict, one per band, each with
        ``name``            preferred model name (must be a ``*_trap`` profile)
        ``t``, ``flux``     the (non-flare) points that were fit
        ``dip_windows``     shared labeled dip windows (same list for all bands)
        ``baseline_params`` the Stage-1 preferred baseline coefficients
        ``t_ref``, ``period``
        ``sigma0``          per-band p2p noise
        ``f``               per-band Stage-1 rescale (sigma_eff = f*sigma0)
        ``init``            (optional) the Stage-1 preferred :class:`ModelResult`,
                            used to assemble the group initial-vs-joint BIC.

    Returns
    -------
    dict with keys ``converged``, ``shared`` (list of per-dip {t0, t0_err, W,
    W_err, T14, T14_err}), ``bands`` (per input band: joint ``dips`` plus its
    joint-stage chi2/BIC at the SAME sigma_eff as Stage-1), and ``group``
    (group-level joint chi2/BIC and, when every ``init`` was supplied, the
    initial chi2/BIC and ``dBIC_tie`` = BIC_joint - BIC_initial).  ``None`` if no
    dip window is populated in any band.

    Per-band joint parameter counts attribute each shared (t0, W) pair as
    ``2 / B_i`` to each of the ``B_i`` bands that cover dip ``i`` (so the
    per-band counts sum to the integer group count).  The group BIC uses the
    rigorous ``K*ln(N_total) - 2*lnL`` with integer K and total points N.
    """
    nbnd = len(bands)
    windows = [tuple(w) for w in bands[0]["dip_windows"]]
    D = len(windows)

    Ts = [np.asarray(b["t"], dtype=float) for b in bands]
    Fs = [np.asarray(b["flux"], dtype=float) for b in bands]
    baselines = [models.baseline_value(b["name"], b["baseline_params"], Ts[i],
                                       b["t_ref"], b["period"])
                 for i, b in enumerate(bands)]
    targets = [baselines[i] - Fs[i] for i in range(nbnd)]   # positive in dip
    sigma_eff = [float(b["f"]) * float(b["sigma0"]) for b in bands]

    masks = [[(Ts[bi] >= w[0]) & (Ts[bi] <= w[1]) for w in windows]
             for bi in range(nbnd)]
    present_bands = [[bi for bi in range(nbnd) if masks[bi][i].any()]
                     for i in range(D)]
    shared_dips = [i for i in range(D) if present_bands[i]]
    if not shared_dips:
        return None
    # --- pack parameters: shared [t0_i, W_i] then per-band [A_{b,i}, r_{b,i}] --
    # The shared half-duration W is capped at k*(dw/2) and the shared mid-time t0
    # is pinned to the labeled window (from the primary band), exactly as in the
    # single-band fit, so a tied dip cannot widen or drift out of its box to
    # merge with a neighbor (see TRAP_W_CAP_FACTOR).
    p0, lo, hi = [], [], []
    shared_idx = {}
    for i in shared_dips:
        x0, x1 = windows[i]
        dw = x1 - x0
        W_cap = max(TRAP_W_CAP_FACTOR * 0.5 * dw, 1e-4)
        t0s, Ws = [], []
        for bi in present_bands[i]:
            m = masks[bi][i]
            t0s.append(float(Ts[bi][m][np.argmin(Fs[bi][m])]))
            Ws.append(max(dw / 2.0, 1e-4))
        shared_idx[i] = (len(p0), len(p0) + 1)
        p0 += [float(np.median(t0s)), float(np.median(Ws))]
        lo += [x0, 1e-5]
        hi += [x1, W_cap]
    band_idx = [dict() for _ in range(nbnd)]
    for bi in range(nbnd):
        for i in shared_dips:
            if bi not in present_bands[i]:
                continue
            m = masks[bi][i]
            local_base = float(np.nanmedian(baselines[bi][m]))
            depth = max(local_base - float(np.nanmin(Fs[bi][m])), 1e-4)
            band_idx[bi][i] = (len(p0), len(p0) + 1)
            p0 += [depth, 0.25]
            lo += [0.0, 1e-4]
            hi += [1.0, 1.0]
    p0 = np.array(p0)
    lo = np.array(lo)
    hi = np.array(hi)

    def resid(p):
        out = []
        for bi in range(nbnd):
            dec = np.zeros_like(Ts[bi])
            for i, idx in band_idx[bi].items():
                t0 = p[shared_idx[i][0]]
                W = p[shared_idx[i][1]]
                A, r = p[idx[0]], p[idx[1]]
                dec = dec + models._trapezoid_dip(Ts[bi], A, t0, W, r * W)
            out.append((dec - targets[bi]) / sigma_eff[bi])
        return np.concatenate(out)

    sol = least_squares(resid, p0, bounds=(lo, hi), method="trf", max_nfev=40000)
    p = sol.x
    converged = bool(sol.success)
    try:
        perr = np.sqrt(np.clip(np.diag(np.linalg.inv(sol.jac.T @ sol.jac)),
                               0.0, np.inf))
    except np.linalg.LinAlgError:
        perr = np.full(p.size, np.nan)

    shared = []
    for i in shared_dips:
        it0, iW = shared_idx[i]
        W, We = float(p[iW]), float(perr[iW])
        shared.append({"dip": i, "t0": float(p[it0]), "t0_err": float(perr[it0]),
                       "W": W, "W_err": We, "T14": 2.0 * W, "T14_err": 2.0 * We})

    band_results = []
    g_chi2 = g_lnL = g_k = g_n = 0.0
    for bi, b in enumerate(bands):
        nb = len(b["baseline_params"])
        n = int(Ts[bi].size)
        dips, k_dip = [], 0.0
        dec = np.zeros_like(Ts[bi])
        for i in [i for i in shared_dips if i in band_idx[bi]]:
            it0, iW = shared_idx[i]
            iA, ir = band_idx[bi][i]
            t0, W, A, r = (float(p[it0]), float(p[iW]),
                           float(p[iA]), float(p[ir]))
            t0e, We, Ae, re = (float(perr[it0]), float(perr[iW]),
                               float(perr[iA]), float(perr[ir]))
            tau = r * W
            tau_e = float(np.hypot(r * We, W * re))
            dec = dec + models._trapezoid_dip(Ts[bi], A, t0, W, tau)
            base0 = float(models.baseline_value(
                b["name"], b["baseline_params"], [t0], b["t_ref"],
                b["period"])[0])
            dur = float(models.trap_fwhm(W, tau))
            dur_e = float(np.hypot(2.0 * We, tau_e))
            dips.append(DipResult(
                A=A, A_err=Ae, t0=t0, t0_err=t0e, w=W, w_err=We,
                delta_pct=float(100.0 * A / base0),
                delta_pct_err=float(100.0 * Ae / base0),
                duration=dur, duration_err=dur_e,
                profile="trap", tau=tau, tau_err=tau_e))
            k_dip += 2.0 + 2.0 / len(present_bands[i])
        resid_b = Fs[bi] - (baselines[bi] - dec)
        chi2_0 = float(np.sum(resid_b ** 2) / b["sigma0"] ** 2)
        f = float(b["f"])
        se = f * float(b["sigma0"])
        chi2 = chi2_0 / f ** 2
        k = nb + k_dip + 1.0
        nu = max(n - k, 1)
        lnL = -0.5 * (chi2 + n * np.log(2.0 * np.pi * se ** 2))
        band_results.append({
            "present_dips": [i for i in shared_dips if i in band_idx[bi]],
            "dips": dips, "n": n, "k": float(k),
            "chi2_red0": float(chi2_0 / nu), "chi2_red": float(chi2 / nu),
            "lnL": float(lnL), "BIC": float(k * np.log(n) - 2.0 * lnL),
            "sigma0": float(b["sigma0"]), "f": f, "sigma_eff": float(se)})
        g_chi2 += chi2
        g_lnL += lnL
        g_k += k
        g_n += n

    group = {"n": int(g_n), "k_joint": float(g_k), "chi2_joint": float(g_chi2),
             "lnL_joint": float(g_lnL),
             "BIC_joint": float(g_k * np.log(g_n) - 2.0 * g_lnL)}
    if all(b.get("init") is not None for b in bands):
        ki = sum(b["init"].k for b in bands)
        li = sum(b["init"].lnL for b in bands)
        ci = sum(b["init"].chi2_red * max(b["init"].n - b["init"].k, 1)
                 for b in bands)
        group.update(
            k_initial=float(ki), lnL_initial=float(li), chi2_initial=float(ci),
            BIC_initial=float(ki * np.log(g_n) - 2.0 * li),
            dBIC_tie=float((g_k * np.log(g_n) - 2.0 * g_lnL)
                           - (ki * np.log(g_n) - 2.0 * li)))

    return {"converged": converged, "shared": shared,
            "bands": band_results, "group": group}
