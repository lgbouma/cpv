import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fit_dips import loaders  # noqa: E402


def test_generic_csv_normalizes_sorts_windows(tmp_path):
    p = tmp_path / "lc.csv"
    # deliberately unsorted, with one point outside the window
    p.write_text("t,f,e\n2.0,20,2\n1.0,10,1\n3.0,30,3\n9.0,90,9\n")
    t, flux, err = loaders.load_lightcurve({
        "loader": "generic_csv", "path": str(p),
        "time_col": "t", "flux_col": "f", "fluxerr_col": "e",
        "window": [0.5, 3.5]})
    assert np.all(np.diff(t) > 0)          # sorted
    assert t.max() <= 3.5 and 9.0 not in t  # windowed
    assert np.isclose(np.median(flux), 1.0)  # normalized to median 1
    assert err is not None


def test_mask_indices_drops_points_and_renormalizes(tmp_path):
    p = tmp_path / "lc.csv"
    p.write_text("t,f\n1.0,10\n2.0,20\n3.0,30\n4.0,40\n5.0,50\n")
    base = {"loader": "generic_csv", "path": str(p),
            "time_col": "t", "flux_col": "f"}
    t0, f0, _ = loaders.load_lightcurve(base)
    t1, f1, _ = loaders.load_lightcurve(dict(base, mask_indices=[0]))
    assert len(t1) == len(t0) - 1
    assert t1[0] == 2.0                      # first (earliest) point dropped
    assert np.isclose(np.median(f1), 1.0)    # re-normalized after masking


def test_fold_period_collapses_cycles_and_resorts(tmp_path):
    p = tmp_path / "lc.csv"
    # three cycles (P=1) of the same two phases (0.2 and 0.7); folding on P=1
    # about fold_ref=0.45 must map every cycle's points onto the base cycle.
    p.write_text("t,f\n0.2,10\n0.7,12\n1.2,10\n1.7,12\n2.2,10\n2.7,12\n")
    base = {"loader": "generic_csv", "path": str(p),
            "time_col": "t", "flux_col": "f"}
    t, flux, _ = loaders.load_lightcurve(
        dict(base, fold_period=1.0, fold_ref=0.45))
    assert np.all(np.diff(t) >= 0)               # re-sorted by folded time
    assert t.min() >= 0.45 - 0.5 - 1e-9          # within +/- P/2 of fold_ref
    assert t.max() <= 0.45 + 0.5 + 1e-9
    # the three cycles collapse to two distinct phases (0.2 and 0.7).
    assert np.allclose(np.unique(np.round(t, 6)), [0.2, 0.7])


def test_unknown_loader_raises():
    try:
        loaders.load_lightcurve({"loader": "nope", "path": "x"})
    except ValueError:
        return
    raise AssertionError("expected ValueError for unknown loader")


def test_loaders_registered():
    for name in ("muscat_csv", "muscat2_fits", "tierras_csv", "keplercam_dat",
                 "fourstar_xls", "tess_fits", "tess_pkl", "generic_csv"):
        assert name in loaders.LOADERS
