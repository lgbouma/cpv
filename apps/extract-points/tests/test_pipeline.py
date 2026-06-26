"""Smoke tests for the digitization pipeline.

Run:  python -m pytest tests/   (or)   python tests/test_pipeline.py
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from extract_points.image_io import load_rgb, dark_mask
from extract_points.panels import detect_flux_boxes
from extract_points.calibrate import calibrate
from extract_points.extract import extract_all, to_dataframe


def test_panel_counts():
    assert len(detect_flux_boxes("f12", load_rgb("f12"))) == 12
    assert len(detect_flux_boxes("f13", load_rgb("f13"))) == 8


def test_calibration_quality():
    for fig in ("f12", "f13"):
        rgb = load_rgb(fig)
        dark = dark_mask(rgb)
        for b in detect_flux_boxes(fig, rgb):
            c = calibrate(dark, b)
            assert len(c.x_ticks_px) >= 3, (fig, b.row, b.col)
            assert c.x_r2 > 0.999, (fig, b.row, b.col, c.x_r2)
            assert c.y_r2 > 0.999, (fig, b.row, b.col, c.y_r2)


def test_extraction_sane():
    # all extraction methods must produce sane, well-sampled curves
    for method in ("per_marker", "cadence", "red_aware"):
        results = extract_all(method)
        assert len(results) == 20, method
        df = to_dataframe(results)
        assert not df.isna().any().any(), method
        # relative flux must sit in the plotted range (I_C ~1.0, IR offset ~0.9)
        assert df.rel_flux.between(0.85, 1.10).all(), method
        # every band has a reasonable number of points (one per ~marker; the
        # sparsest are partial-coverage non-CSV panels)
        for r in results:
            for bd in r.bands:
                assert bd.bjd.size >= 15, (method, r.epoch, bd.band, bd.bjd.size)
                # points are time-ordered and span a non-trivial range
                assert np.all(np.diff(bd.bjd) >= 0)
                assert np.ptp(bd.bjd) > 0.05


def test_csv_epochs_present():
    for method in ("per_marker", "cadence", "red_aware"):
        results = extract_all(method)
        csv_epochs = {r.epoch for r in results if r.in_csv}
        assert len(csv_epochs) == 9, method


if __name__ == "__main__":
    test_panel_counts()
    test_calibration_quality()
    test_extraction_sane()
    test_csv_epochs_present()
    print("All smoke tests passed.")
