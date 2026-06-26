"""Tests for complexrotators.fit_dips_beta.load_depths category filtering."""
from os.path import dirname, join

import pandas as pd

from complexrotators import fit_dips_beta as fdb

CSV_PATH = join(
    dirname(dirname(__file__)), "papers", "Bouma_2026_cgcd",
    "tables", "dipdepths_thiswork_and_literature.csv",
)


def test_combined_is_unchanged_by_none_categories():
    df = fdb.load_depths(CSV_PATH, categories=None)
    assert len(df) == 57
    assert set(df["category"]) == {"This work", "Literature"}
    assert df["star"].nunique() == 6


def test_thiswork_only():
    df = fdb.load_depths(CSV_PATH, categories=["This work"])
    assert len(df) == 27
    assert set(df["category"]) == {"This work"}
    assert set(df["star"]) == {
        "LP 12-502", "TIC 262400835", "TIC 300651846"
    }


def test_literature_only():
    df = fdb.load_depths(CSV_PATH, categories=["Literature"])
    assert len(df) == 30
    assert set(df["category"]) == {"Literature"}
    assert set(df["star"]) == {
        "PTFO 8-8695", "TIC 20178925", "TIC 435899024"
    }


def test_row_ids_anchored_to_full_table():
    """A this-work row keeps the row_id it has in the combined table."""
    full = fdb.load_depths(CSV_PATH, categories=None)
    tw = fdb.load_depths(CSV_PATH, categories=["This work"])
    full_tw_ids = set(full.loc[full["category"] == "This work", "row_id"])
    assert set(tw["row_id"]) == full_tw_ids


def test_categories_and_exclusions_compose():
    """Excluding a this-work row_id still works under a category filter."""
    tw = fdb.load_depths(CSV_PATH, categories=["This work"])
    drop = int(tw["row_id"].iloc[0])
    tw2 = fdb.load_depths(
        CSV_PATH, exclude_row_ids=[drop], categories=["This work"]
    )
    assert drop not in set(tw2["row_id"])
    assert len(tw2) == len(tw) - 1
