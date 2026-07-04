"""Shared persistence for coincidence labels.

A *label* marks an x-location in one mosaic column where a feature is judged to
coincide in BOTH the Hα-EW (top) and Δflux (bottom) panels:

  kind='dip'    local minima coincide  -> drawn as a '|' tick
  kind='flare'  local maxima coincide  -> drawn as a 'v' tick

Labels are stored in a single long-form CSV (columns: figure, colkey, kind, x)
so the interactive app writes them and the two plotting drivers read them.

x is stored in each column's native x-axis units: hours-from-start for
mode='time' columns, orbital phase for mode='phase' columns.

This module deliberately depends only on the standard library + pandas so it can
be imported cheaply from both the app and the drivers.
"""
from __future__ import annotations

import os
from os.path import dirname, abspath, join

import pandas as pd

# output/coincidence_labels.csv, next to this package
APPDIR = dirname(dirname(abspath(__file__)))
LABELS_PATH = join(APPDIR, "output", "coincidence_labels.csv")

KINDS = ("dip", "flare")
_COLUMNS = ["figure", "colkey", "kind", "x"]


def colkey(col: dict) -> str:
    """Stable identifier for a mosaic column dict (as defined in the drivers).

    Content-based (targetid + date[s]) so it survives column reordering.
    """
    tid = col["targetid"]
    if col["mode"] == "time":
        return f"{tid}_{col['utcdatestr']}"
    return f"{tid}_" + "+".join(col["utcdatestrs"])


def load_labels() -> pd.DataFrame:
    """Return all labels as a DataFrame (empty with the right columns if none)."""
    if not os.path.exists(LABELS_PATH):
        return pd.DataFrame(columns=_COLUMNS)
    df = pd.read_csv(LABELS_PATH)
    for c in _COLUMNS:
        if c not in df.columns:
            df[c] = pd.Series(dtype="float" if c == "x" else "object")
    return df[_COLUMNS]


def save_labels(df: pd.DataFrame) -> None:
    """Write labels back to LABELS_PATH, sorted for a stable diff."""
    os.makedirs(dirname(LABELS_PATH), exist_ok=True)
    df = df[_COLUMNS].sort_values(["figure", "colkey", "kind", "x"])
    df.to_csv(LABELS_PATH, index=False)


def get_column_labels(df: pd.DataFrame, figure: str, ck: str, kind: str):
    """Sorted list of x-values for one (figure, colkey, kind)."""
    if df.empty:
        return []
    sel = (df.figure == figure) & (df.colkey == ck) & (df.kind == kind)
    return sorted(df.loc[sel, "x"].astype(float).tolist())
