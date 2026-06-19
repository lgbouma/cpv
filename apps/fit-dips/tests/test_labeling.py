import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fit_dips.labeling import LabelState  # noqa: E402


def test_add_and_masks():
    t = np.linspace(0, 1, 101)
    s = LabelState()
    s.add_dip(0.2, 0.3)
    s.add_flare(0.5, 0.55)
    in_dip, in_flare, oot = s.masks(t)
    assert in_dip.sum() > 0 and in_flare.sum() > 0
    # flare and dip are disjoint; everything sums to all points
    assert np.all(in_dip + in_flare + oot == 1)
    assert not np.any(in_dip & in_flare)


def test_flare_overrides_dip_overlap():
    t = np.linspace(0, 1, 101)
    s = LabelState()
    s.add_dip(0.2, 0.6)
    s.add_flare(0.3, 0.4)
    in_dip, in_flare, oot = s.masks(t)
    # points inside the flare window are NOT counted as in-dip
    overlap = (t >= 0.3) & (t <= 0.4)
    assert not np.any(in_dip & overlap)
    assert np.all(in_flare == overlap)


def test_undo_and_reset():
    s = LabelState()
    s.add_dip(0.1, 0.2)
    s.add_flare(0.3, 0.4)
    s.undo()
    assert len(s.flares) == 0 and len(s.dips) == 1
    s.reset()
    assert len(s.dips) == 0 and len(s.flares) == 0 and s.order == []
