"""
Dip + baseline models for CPV dimming events.

The combined model evaluated over the non-flare points is

    F(t) = B(t; theta_b) - sum_i  D_i(t; profile)

where B is one of the baseline families (polynomial or Fourier) and D_i is a dip
profile.  Two dip profiles are supported and compared by BIC:

  "sech" : A * sech( (t - t0) / w )                       (Yu+2015 Eq. 1; 3 par)
  "trap" : a flat-bottomed transit (symmetric trapezoid)  (4 par)
           full-depth A for |t-t0| <= W - tau; linear ingress/egress to 0 at
           |t-t0| = W.  Nests a box (tau -> 0) and a triangle (tau -> W).

A *model name* encodes both the baseline family and the dip profile: the sech
profile uses the bare baseline name ("poly2"), the trapezoid appends a suffix
("poly2_trap").  parse_model_name() / make_model_name() convert between the two.

Conventions
-----------
- Reported depth   : delta_i = A_i / B(t0_i), in percent (A_i is the full-depth
  decrement at mid-dip for BOTH profiles, so depths are directly comparable).
- Reported duration: FWHM (width at half depth).
    sech       : 2*arccosh(2)*w ~= 2.6339*w.
    trapezoid  : 2*W - tau  (half-depth falls at the ingress/egress midpoint).
- Reported midtime : t0_i.
"""
import numpy as np

# FWHM of sech(x): sech(x)=1/2 at x = arccosh(2); full width = 2*arccosh(2).
ARCCOSH2 = np.arccosh(2.0)
FWHM_OVER_W = 2.0 * ARCCOSH2  # ~= 2.6339

# All baseline families that exist.  The DEFAULT pool (used by fit_all_models)
# is the lower-order subset -- poly3/poly4/fourier3/fourier4 tend to overfit the
# short single-night windows -- but every family remains available by name
# (e.g. via fit_one_model) for diagnostics.
BASELINE_FAMILIES = [
    "poly1", "poly2", "poly3", "poly4",
    "fourier1", "fourier2", "fourier3", "fourier4",
]
DEFAULT_BASELINE_FAMILIES = ["poly1", "poly2", "fourier1", "fourier2"]

# Dip profiles that exist.  The DEFAULT pool (used by fit_all_models) is the
# flat-bottomed transit ("trap") ONLY: the cusped sech systematically misfit and
# inflated the depths of these dips, so the default model comparison is now just
# the trapezoid against each local-trend baseline.  The sech remains available by
# name (e.g. fit_one_model("poly1", ...)) for diagnostics / comparison.
DIP_PROFILES = ["sech", "trap"]
DEFAULT_DIP_PROFILES = ["trap"]
_TRAP_SUFFIX = "_trap"


def n_dip_params(profile):
    """Number of free parameters per dip for a profile."""
    return 4 if profile == "trap" else 3


def parse_model_name(name):
    """Split a model name into (baseline_family, dip_profile).

    "poly2"      -> ("poly2", "sech")
    "poly2_trap" -> ("poly2", "trap")
    """
    if name.endswith(_TRAP_SUFFIX):
        return name[: -len(_TRAP_SUFFIX)], "trap"
    return name, "sech"


def make_model_name(baseline, profile):
    """Compose a model name from a baseline family and a dip profile."""
    return baseline + _TRAP_SUFFIX if profile == "trap" else baseline


def default_model_names():
    """Default model pool: each DEFAULT baseline x each DEFAULT dip profile.

    With DEFAULT_DIP_PROFILES = ["trap"] this is the trapezoid combined with each
    local-trend baseline (poly1, poly2, fourier1, fourier2).
    """
    return [make_model_name(b, p)
            for b in DEFAULT_BASELINE_FAMILIES for p in DEFAULT_DIP_PROFILES]


def sech(x):
    """Numerically stable hyperbolic secant, 2/(e^x + e^-x) = 1/cosh(x)."""
    return 1.0 / np.cosh(np.asarray(x, dtype=float))


def _sech_dip(t, A, t0, w):
    return A * sech((t - t0) / w)


def _trapezoid_dip(t, A, t0, W, tau):
    """Flat-bottomed (symmetric trapezoid) transit decrement, >= 0.

    Full depth ``A`` for |t-t0| <= W - tau; linear ingress/egress down to 0 at
    |t-t0| = W (so W is the half total duration, tau the ingress/egress
    duration).  tau is clipped to (0, W] for safety.
    """
    t = np.asarray(t, dtype=float)
    W = float(W)
    tau = float(np.clip(tau, 1e-12, W if W > 0 else 1e-12))
    x = np.abs(t - t0)
    flat = W - tau
    ramp = A * (W - x) / tau           # full-depth at x=flat, 0 at x=W
    return np.where(x <= flat, A, np.where(x <= W, ramp, 0.0))


def dip_model(t, dip_params, profile="sech"):
    """Sum of dip profiles.

    Parameters
    ----------
    t : array
    dip_params : flat array, ``n_dip_params(profile)`` values per dip:
        sech -> [A, t0, w, ...]; trap -> [A, t0, W, tau, ...].
    profile : "sech" or "trap"

    Returns
    -------
    array of the (positive) flux *decrement* contributed by all dips.
    """
    t = np.asarray(t, dtype=float)
    out = np.zeros_like(t)
    npp = n_dip_params(profile)
    p = np.asarray(dip_params, dtype=float).reshape(-1, npp)
    for row in p:
        if profile == "trap":
            A, t0, W, tau = row
            out = out + _trapezoid_dip(t, A, t0, W, tau)
        else:
            A, t0, w = row
            out = out + _sech_dip(t, A, t0, w)
    return out


def n_baseline_params(name):
    """Number of free baseline parameters for a (possibly profile-suffixed) name."""
    name, _ = parse_model_name(name)
    if name.startswith("poly"):
        return int(name[len("poly"):]) + 1          # poly1 -> 2 (c0, c1)
    if name.startswith("fourier"):
        return 2 * int(name[len("fourier"):]) + 1    # fourier1 -> 3 (c0, a1, b1)
    raise ValueError(f"unknown baseline family: {name}")


def baseline_design(name, t, t_ref, period=None):
    """Design matrix X such that baseline = X @ params.

    Time is centered on ``t_ref`` for numerical conditioning.  Fourier families
    require ``period`` (days).  A trailing dip-profile suffix is ignored.
    """
    name, _ = parse_model_name(name)
    t = np.asarray(t, dtype=float)
    tau = t - t_ref
    if name.startswith("poly"):
        deg = int(name[len("poly"):])
        return np.vstack([tau ** k for k in range(deg + 1)]).T
    if name.startswith("fourier"):
        N = int(name[len("fourier"):])
        if period is None:
            raise ValueError(f"{name} requires a period")
        cols = [np.ones_like(t)]
        ph = 2.0 * np.pi * tau / period
        for k in range(1, N + 1):
            cols.append(np.cos(k * ph))
            cols.append(np.sin(k * ph))
        return np.vstack(cols).T
    raise ValueError(f"unknown baseline family: {name}")


def baseline_value(name, params, t, t_ref, period=None):
    """Evaluate the baseline at times ``t``."""
    X = baseline_design(name, t, t_ref, period)
    return X @ np.asarray(params, dtype=float)


def full_model(name, baseline_params, dip_params, t, t_ref, period=None,
               profile=None):
    """Baseline minus dips.  Dip profile is taken from ``name`` unless given."""
    if profile is None:
        _, profile = parse_model_name(name)
    return (
        baseline_value(name, baseline_params, t, t_ref, period)
        - dip_model(t, dip_params, profile)
    )


def w_to_fwhm(w):
    """Convert sech width parameter to FWHM duration."""
    return FWHM_OVER_W * np.asarray(w, dtype=float)


def trap_fwhm(W, tau):
    """FWHM (width at half depth) of the trapezoid dip = 2*W - tau."""
    return 2.0 * np.asarray(W, dtype=float) - np.asarray(tau, dtype=float)
