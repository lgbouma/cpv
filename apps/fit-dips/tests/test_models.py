import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fit_dips import models  # noqa: E402


def test_sech_peak_and_symmetry():
    assert np.isclose(models.sech(0.0), 1.0)
    assert np.isclose(models.sech(1.3), models.sech(-1.3))


def test_fwhm_relation():
    # sech(x) = 1/2 at x = arccosh(2); full width = 2*arccosh(2)*w.
    w = 0.01
    fwhm = models.w_to_fwhm(w)
    half_x = fwhm / 2.0 / w
    assert np.isclose(models.sech(half_x), 0.5)
    assert np.isclose(fwhm, 2.6339 * w, atol=1e-3)


def test_baseline_param_counts():
    assert models.n_baseline_params("poly1") == 2
    assert models.n_baseline_params("poly4") == 5
    assert models.n_baseline_params("fourier1") == 3
    assert models.n_baseline_params("fourier4") == 9


def test_baseline_design_shapes():
    t = np.linspace(0, 1, 50)
    assert models.baseline_design("poly3", t, 0.5).shape == (50, 4)
    assert models.baseline_design("fourier2", t, 0.5, period=0.3).shape == (50, 5)


def test_dip_model_depth_at_center():
    t = np.array([5.0])
    val = models.dip_model(t, [0.02, 5.0, 0.01])
    assert np.isclose(val[0], 0.02)


def test_full_model_subtracts_dip():
    t = np.linspace(0, 1, 11)
    base = models.baseline_value("poly1", [1.0, 0.0], t, 0.5)
    full = models.full_model("poly1", [1.0, 0.0], [0.1, 0.5, 0.05], t, 0.5)
    assert np.all(full <= base + 1e-12)


def test_parse_and_make_model_name():
    assert models.parse_model_name("poly2") == ("poly2", "sech")
    assert models.parse_model_name("poly2_trap") == ("poly2", "trap")
    assert models.parse_model_name("fourier1_trap") == ("fourier1", "trap")
    assert models.make_model_name("poly2", "sech") == "poly2"
    assert models.make_model_name("poly2", "trap") == "poly2_trap"
    # baseline helpers ignore the profile suffix
    assert models.n_baseline_params("poly2_trap") == 3
    assert models.baseline_design("poly2_trap", np.linspace(0, 1, 5),
                                  0.5).shape == (5, 3)


def test_n_dip_params_and_default_pool():
    assert models.n_dip_params("sech") == 3
    assert models.n_dip_params("trap") == 4
    # default pool drops poly3/4 + fourier3/4 AND the sech profile -> trap only.
    assert models.DEFAULT_DIP_PROFILES == ["trap"]
    assert models.default_model_names() == [
        "poly1_trap", "poly2_trap", "fourier1_trap", "fourier2_trap"]


def test_trapezoid_flat_bottom_ramp_and_zero():
    # A=0.05, t0=0, W=0.04 (half total), tau=0.01 ingress => flat |x|<=0.03.
    t = np.array([0.0, 0.025, 0.035, 0.05, 0.06])
    d = models.dip_model(t, [0.05, 0.0, 0.04, 0.01], "trap")
    assert np.isclose(d[0], 0.05)            # center: full depth
    assert np.isclose(d[1], 0.05)            # still on the flat bottom
    assert np.isclose(d[2], 0.025)           # ingress midpoint -> half depth
    assert np.isclose(d[3], 0.0)             # first/fourth contact
    assert np.isclose(d[4], 0.0)             # outside
    assert np.all(d >= 0.0)
    # FWHM (width at half depth) = 2W - tau.
    assert np.isclose(models.trap_fwhm(0.04, 0.01), 0.07)


def test_trapezoid_box_and_triangle_limits():
    t = np.linspace(-0.05, 0.05, 2001)
    box = models.dip_model(t, [0.03, 0.0, 0.02, 1e-9], "trap")     # tau->0
    assert np.isclose(box[np.argmin(np.abs(t - 0.019))], 0.03)     # flat to edge
    tri = models.dip_model(t, [0.03, 0.0, 0.02, 0.02], "trap")     # tau->W
    assert np.isclose(tri[np.argmin(np.abs(t))], 0.03)             # peak only
    assert tri[np.argmin(np.abs(t - 0.01))] < 0.03                 # already ramping
