#!/usr/bin/env python
"""Smoke tests for label-coincidences.

Run:  python tests/test_labels.py   (with the cpv env)
"""
import os
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from label_coincidences import labels as L
from label_coincidences.data import figure_columns, build_column_series


def test_colkey():
    time_col = dict(mode="time", targetid="LP_12-502", utcdatestr="20241104")
    phase_col = dict(mode="phase", targetid="LP_12-502",
                     utcdatestrs=["20231111", "20231112"])
    assert L.colkey(time_col) == "LP_12-502_20241104"
    assert L.colkey(phase_col) == "LP_12-502_20231111+20231112"
    print("ok  colkey")


def test_roundtrip(tmp_path):
    orig = L.LABELS_PATH
    L.LABELS_PATH = os.path.join(tmp_path, "coincidence_labels.csv")
    try:
        import pandas as pd
        assert L.load_labels().empty  # missing file -> empty
        df = pd.DataFrame(
            [("mosaic", "DG_CVn_20240428", "dip", 1.5),
             ("mosaic", "DG_CVn_20240428", "flare", 3.0)],
            columns=["figure", "colkey", "kind", "x"],
        )
        L.save_labels(df)
        got = L.load_labels()
        assert len(got) == 2
        assert L.get_column_labels(got, "mosaic", "DG_CVn_20240428", "dip") == [1.5]
        assert L.get_column_labels(got, "mosaic", "DG_CVn_20240428", "flare") == [3.0]
        assert L.get_column_labels(got, "mosaic", "DG_CVn_20240428", "dip") != [3.0]
        print("ok  save/load roundtrip")
    finally:
        L.LABELS_PATH = orig


def test_data_and_session():
    n = 0
    for fig in ("mosaic", "supp"):
        for col in figure_columns(fig):
            s = build_column_series(fig, col)
            assert len(s["ew"]["x"]) > 0
            assert len(s["flux"]["x"]) > 0
            assert len(s["snap_x"]) == len(s["ew"]["x"])
            assert s["xlim"][0] < s["xlim"][1]
            n += 1
    print(f"ok  build_column_series for {n} columns")

    from label_coincidences.interactive import LabelSession
    sess = LabelSession("supp")
    sess.state[0]["dip"].append(sess._snap(0, 2.0))
    sess.state[0]["flare"].append(sess._snap(0, 5.0))
    sess._redraw_col(0)
    # both kinds now present -> one artist per kind per panel = 4
    assert len(sess.marker_artists[0]) == 4
    print("ok  LabelSession builds + draws markers")


if __name__ == "__main__":
    test_colkey()
    with tempfile.TemporaryDirectory() as d:
        test_roundtrip(d)
    test_data_and_session()
    print("\nAll tests passed.")
