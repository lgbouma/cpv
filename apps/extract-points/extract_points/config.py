"""Per-panel metadata, hand-read from the two figures.

For each flux panel we record the epoch label, the BJD offset shown in the
x-axis label, the photometric bands present (optical I_C plotted as filled
circles, IR plotted as crosses), and `x_first`: the value of the leftmost
major tick on the bottom (BJD - offset) axis. Major ticks are spaced by
`X_STEP` = 0.1; the number of ticks is recovered automatically by tick
detection, so only the leftmost value is needed here.

The y-axis (relative flux) has major ticks at 1.0, 0.9, ... (top tick = 1.0,
step -0.1); this is uniform across all panels and handled in calibrate.py.
"""
from __future__ import annotations

X_STEP = 0.1
Y_TOP = 1.0
Y_STEP = 0.1

# key: (figure, row, col)
PANELS: dict[tuple[str, int, int], dict] = {
    # ----- f12 (4 rows x 3 cols) -----
    ("f12", 0, 0): dict(epoch="2014 Feb 23", offset=2456711, bands=["I_C", "K_s"], x_first=0.9),
    ("f12", 0, 1): dict(epoch="2014 Dec 27", offset=2457019, bands=["I_C", "K_s"], x_first=0.0),
    ("f12", 0, 2): dict(epoch="2015 Jan 10", offset=2457032, bands=["I_C", "K_s"], x_first=0.9),
    ("f12", 1, 0): dict(epoch="2015 Jan 23", offset=2457045, bands=["I_C", "J"], x_first=0.9),
    ("f12", 1, 1): dict(epoch="2015 Feb 10", offset=2457063, bands=["I_C", "H"], x_first=0.9),
    ("f12", 1, 2): dict(epoch="2015 Feb 14", offset=2457067, bands=["I_C", "K_s"], x_first=0.9),
    ("f12", 2, 0): dict(epoch="2015 Feb 23", offset=2457076, bands=["I_C", "K_s"], x_first=0.9),
    ("f12", 2, 1): dict(epoch="2015 Oct 17", offset=2457313, bands=["I_C"], x_first=0.1),
    ("f12", 2, 2): dict(epoch="2016 Nov 3", offset=2457696, bands=["I_C", "K_s"], x_first=0.1),
    ("f12", 3, 0): dict(epoch="2016 Nov 25", offset=2457718, bands=["I_C", "J"], x_first=0.1),
    ("f12", 3, 1): dict(epoch="2016 Nov 29", offset=2457721, bands=["I_C", "J"], x_first=1.0),
    ("f12", 3, 2): dict(epoch="2016 Dec 1", offset=2457724, bands=["I_C", "J"], x_first=0.1),
    # ----- f13 (rows of 3,3,2) -----
    ("f13", 0, 0): dict(epoch="2016 Dec 9", offset=2457732, bands=["I_C", "J"], x_first=0.0),
    ("f13", 0, 1): dict(epoch="2017 Oct 3", offset=2458030, bands=["I_C", "J"], x_first=0.1),
    ("f13", 0, 2): dict(epoch="2018 Feb 8", offset=2458158, bands=["I_C", "J"], x_first=-0.1),
    ("f13", 1, 0): dict(epoch="2018 Nov 7", offset=2458430, bands=["I_C"], x_first=0.1),
    ("f13", 1, 1): dict(epoch="2018 Nov 9", offset=2458432, bands=["I_C"], x_first=0.1),
    ("f13", 1, 2): dict(epoch="2018 Nov 10", offset=2458433, bands=["I_C"], x_first=0.1),
    ("f13", 2, 0): dict(epoch="2018 Dec 29", offset=2458482, bands=["I_C", "J"], x_first=-0.1),
    ("f13", 2, 1): dict(epoch="2018 Dec 30", offset=2458482, bands=["I_C", "J"], x_first=0.9),
}

# Epochs that overlap dipdepths_thiswork_and_literature.csv (ref=6, PTFO 8-8695)
CSV_EPOCHS = {
    "2014 Feb 23", "2014 Dec 27", "2015 Jan 10", "2015 Jan 23", "2015 Feb 10",
    "2015 Feb 14", "2015 Feb 23", "2016 Nov 25", "2016 Dec 9",
}
