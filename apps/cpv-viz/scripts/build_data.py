"""Build cpv-viz data assets from the literature compilation table."""

from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable

import numpy as np
from astropy.time import Time

ROOT_DIR = Path(__file__).resolve().parents[3]
SOURCE_CSV = ROOT_DIR / (
    "results/tables/jan2026_compilation/"
    "20240304_CPV_lit_compilation_R16_S17_S18_B20_S21_Z19_G22_P23_B24_TIC8_obs_supplemented.csv"
)
OUTPUT_DIR = ROOT_DIR / "apps/cpv-viz/data"
BASE_NAME = "cpv_lit_subset"
MANIFEST_JSON = OUTPUT_DIR / "cpv_lit_subset_manifest.json"
LATEST_JSON = OUTPUT_DIR / f"{BASE_NAME}.json"

COLUMNS = [
    "original_id",
    "ticid",
    "sectors",
    "N_sectors",
    "distance_pc",
    "TESSMAG",
    "cluster",
    "bibcode",
    "period_hr",
    "quality",
    "telescope",
]

FLOAT_COLUMNS = {"distance_pc", "TESSMAG", "period_hr"}
INT_COLUMNS = {"ticid", "quality", "N_sectors"}


def _parse_float(value: str) -> float | None:
    """Parse a float value from a CSV field.

    Args:
        value: Raw string from the CSV.

    Returns:
        Parsed float, or None if the value is empty or NaN-like.
    """
    if value is None:
        return None
    cleaned = value.strip()
    if cleaned == "" or cleaned.lower() == "nan":
        return None
    try:
        return float(cleaned)
    except ValueError:
        return None


def _parse_int(value: str) -> int | None:
    """Parse an integer value from a CSV field.

    Args:
        value: Raw string from the CSV.

    Returns:
        Parsed integer, or None if the value is empty or invalid.
    """
    if value is None:
        return None
    cleaned = value.strip()
    if cleaned == "" or cleaned.lower() == "nan":
        return None
    try:
        return int(float(cleaned))
    except ValueError:
        return None


def _parse_value(column: str, value: str) -> Any:
    """Parse a CSV value based on its column type.

    Args:
        column: Column name.
        value: Raw string from the CSV.

    Returns:
        Parsed value with the appropriate type.
    """
    if column in FLOAT_COLUMNS:
        return _parse_float(value)
    if column in INT_COLUMNS:
        return _parse_int(value)
    if value is None:
        return None
    cleaned = value.strip()
    if cleaned == "" or cleaned.lower() == "nan":
        return None
    return cleaned


def _format_display_date(date_value: datetime) -> str:
    """Format a date for the UI subtitle.

    Args:
        date_value: Datetime to format.

    Returns:
        Human-friendly date string (e.g., "Feb 2, 2026").
    """
    month_names = [
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ]
    month_name = month_names[date_value.month - 1]
    return f"{month_name} {date_value.day}, {date_value.year}"


def _current_sector_number() -> int | None:
    """Return the current TESS sector number based on midtimes.

    Returns:
        The latest sector whose midpoint is not in the future, or None if
        unavailable.
    """
    try:
        from tess_stars2px import TESS_Spacecraft_Pointing_Data as tspd
    except ImportError:
        return None

    sectors = np.array(tspd.sectors, dtype=int)
    times = Time(tspd.midtimes, format="jd")
    now = Time.now()

    in_past = times <= now
    if not np.any(in_past):
        return None

    return int(np.max(sectors[in_past]))


def _load_rows(rows: Iterable[Dict[str, str]]) -> list[Dict[str, Any]]:
    """Extract and parse the columns needed for cpv-viz.

    Args:
        rows: Iterable of CSV rows.

    Returns:
        List of parsed row dictionaries.
    """
    output_rows: list[Dict[str, Any]] = []
    for row in rows:
        parsed_row = {column: _parse_value(column, row.get(column)) for column in COLUMNS}
        output_rows.append(parsed_row)
    return output_rows


def build_data() -> Path:
    """Build the cpv-viz JSON data file.

    Returns:
        Path to the generated JSON file.

    Raises:
        FileNotFoundError: If the source CSV cannot be found.
    """
    if not SOURCE_CSV.exists():
        raise FileNotFoundError(f"Source CSV not found: {SOURCE_CSV}")

    with SOURCE_CSV.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = _load_rows(reader)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    generated_at = datetime.now()
    date_stamp = generated_at.strftime("%Y%m%d")
    output_json = OUTPUT_DIR / f"{BASE_NAME}_{date_stamp}.json"

    with output_json.open("w", encoding="utf-8") as handle:
        json.dump(rows, handle, indent=2)

    with LATEST_JSON.open("w", encoding="utf-8") as handle:
        json.dump(rows, handle, indent=2)

    current_sector = _current_sector_number()
    manifest = {
        "generated_date": generated_at.strftime("%Y-%m-%d"),
        "display_date": _format_display_date(generated_at),
        "data_file": output_json.name,
        "current_sector": current_sector,
    }

    with MANIFEST_JSON.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    return output_json


def main() -> None:
    """Run the data build and print the output path."""
    output_path = build_data()
    print(f"Wrote {output_path}")
    print(f"Wrote {MANIFEST_JSON}")
    print(f"Wrote {LATEST_JSON}")


if __name__ == "__main__":
    main()
