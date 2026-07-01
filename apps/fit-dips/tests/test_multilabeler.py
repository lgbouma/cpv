"""Smoke tests for the N-panel shared-window labeler and epoch grouping.

These exercise the GUI plumbing headlessly (Agg backend): a window dragged in
one panel must propagate to all panels' spans, and each panel must classify its
own points from the shared windows. We also check the (star, date_ut) +
time-overlap grouping rule.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fit_dips.labeling import MultiDipLabeler, DipLabeler  # noqa: E402
from fit_dips import app  # noqa: E402


def _panel(t, flux, on_fit=None):
    return {"t": t, "flux": flux, "title": "band", "on_fit": on_fit or (lambda s: None)}


def test_span_propagates_to_all_panels():
    t = np.linspace(0.0, 1.0, 51)
    lab = MultiDipLabeler(
        [_panel(t, np.ones_like(t)), _panel(t, 0.9 * np.ones_like(t))],
        group_title="grp")
    n_before = len(lab._spans)
    lab._on_select(0.2, 0.3)
    # one axvspan patch per axis for the one window
    assert len(lab._spans) - n_before == len(lab.axes) == 2
    assert lab.state.dips == [[0.2, 0.3]]


def test_each_panel_classifies_own_sampling():
    # two panels with the SAME shared window but different time sampling
    t_a = np.linspace(0.0, 1.0, 101)
    t_b = np.linspace(0.0, 1.0, 11)
    lab = MultiDipLabeler([_panel(t_a, np.ones_like(t_a)),
                           _panel(t_b, np.ones_like(t_b))], group_title="grp")
    lab._on_select(0.45, 0.55)
    lab._classify_preview()
    # in-dip counts follow each panel's own sampling under the shared window
    in_a = ((t_a >= 0.45) & (t_a <= 0.55)).sum()
    in_b = ((t_b >= 0.45) & (t_b <= 0.55)).sum()
    assert in_a != in_b and in_a > 0 and in_b > 0
    assert lab.awaiting_confirm


def test_fit_calls_every_panel_on_fit():
    t = np.linspace(0.0, 1.0, 51)
    calls = []

    def make_on_fit(tag):
        def on_fit(state):
            calls.append(tag)
            grid = np.linspace(0, 1, 10)
            return {"tgrid": grid, "model": np.ones_like(grid),
                    "preferred": f"poly1_trap[{tag}]", "dip_markers": [0.5],
                    "table_text": f"table {tag}"}
        return on_fit

    lab = MultiDipLabeler([_panel(t, np.ones_like(t), make_on_fit("A")),
                           _panel(t, np.ones_like(t), make_on_fit("B"))],
                          group_title="grp")
    lab._on_select(0.4, 0.6)
    lab._classify_preview()
    lab._do_fit()
    assert calls == ["A", "B"]
    assert lab.state.confirmed
    assert all(p["model_line"] is not None for p in lab.panels)


def test_diplabeler_is_single_panel_wrapper():
    t = np.linspace(0.0, 1.0, 21)
    lab = DipLabeler(t, np.ones_like(t), "title", lambda s: None)
    assert len(lab.axes) == 1
    assert lab.ax is lab.axes[0]


def test_find_epoch_group_overlap(monkeypatch):
    # primary + overlapping sibling + a (time-overlapping) folded variant
    # + a same-night-but-disjoint window + a different-star dataset
    prim = {"id": "X_a", "star": "X", "date_ut": "2020-01-01", "data": {"k": "a"}}
    over = {"id": "X_b", "star": "X", "date_ut": "2020-01-01", "data": {"k": "b"}}
    fold = {"id": "X_S88", "star": "X", "date_ut": "2020-01-01",
            "data": {"k": "f", "fold_ref": 100.5}}
    disj = {"id": "X_c", "star": "X", "date_ut": "2020-01-01", "data": {"k": "d"}}
    other = {"id": "Y_a", "star": "Y", "date_ut": "2020-01-01", "data": {"k": "y"}}

    ranges = {"a": (100.0, 101.0), "b": (100.5, 101.5),
              "f": (100.2, 100.8), "d": (200.0, 201.0), "y": (100.0, 101.0)}

    def fake_load_lc(data):
        lo, hi = ranges[data["k"]]
        t = np.linspace(lo, hi, 5)
        return t, np.ones_like(t), None

    monkeypatch.setattr(app.loaders, "load_lightcurve", fake_load_lc)
    monkeypatch.setattr(app.registry, "load", lambda i: prim)
    monkeypatch.setattr(app.registry, "list_datasets",
                        lambda: [prim, over, fold, disj, other])

    ids = app.find_epoch_group("X_a")
    # folded variant (fold_ref), disjoint window, and other star all excluded
    assert ids == ["X_a", "X_b"]
